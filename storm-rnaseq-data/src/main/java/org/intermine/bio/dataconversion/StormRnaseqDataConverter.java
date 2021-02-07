package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 * 
 * @author
 */
public class StormRnaseqDataConverter extends BioDirectoryConverter
{
    private static final String DATASET_TITLE = "STORM RNA-Seq Data";
    private static final String DATA_SOURCE_NAME = "STORM RNA-Seq Data";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private Map<String, String> genes = new HashMap<String, String>();
    private ArrayList<String> conditions = new ArrayList<String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormRnaseqDataConverter.class);

    private String organismIdentifier;
    public StormRnaseqDataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    public void process(File dataDir) throws Exception {
        organismIdentifier = getOrganism(TAXON_ID);

        Map<String, File> directories = readDirectoriesInDir(dataDir);

        for (Map.Entry<String, File> entry : directories.entrySet()) {
            File theFolder = entry.getValue();
            String experimentFolderName = entry.getKey();
            Map<String, File> files = readFilesInDir(theFolder);
            for (Map.Entry<String, File> fileEntry : files.entrySet()) {
                String cellLineFolderName = fileEntry.getKey();
                String experimentName = experimentFolderName + "-" + cellLineFolderName;

                String rootFolder = theFolder.getAbsolutePath();

                // Sample info
                String sampleInfoPath = rootFolder.concat(cellLineFolderName + "_sample_info.csv")
                processRNASeqExperimentSampleInfo(new FileReader(sampleInfoPath), experimentName);

                // Gene counts
                String countsPath = rootFolder.concat("merged_gene_counts.txt")
                processRNASeqExperimentGeneCount(new FileReader(countsPath), experimentName);

                break;
            }
        }

    }

    private void processRNASeqExperimentSampleInfo(Reader reader, String experimentName) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        String[] firstLine = (String[]) lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String run = line[0];
            String sample = line[1];
            String condition = line[2];
            String cellLine = line[3];
            String IFN_gamma = line[4];
            String compound = line[5];
            String concentration = line[6];
            String replicate = line[7];
            String timepoint = line[8];

            Item IntegratedItem = createItem("RNASeqExperimentSampleInfo");

            if(!StringUtils.isEmpty(run)) {
                IntegratedItem.setAttribute("run", run);
            } else {
                continue;
            }

            if(!StringUtils.isEmpty(sample)) {
                IntegratedItem.setAttribute("sample", sample);
            }

            if(!StringUtils.isEmpty(condition)) {
                IntegratedItem.setAttribute("condition", condition);

                if(!condition.contains("DMSO")) {
                    if(!conditions.contains(condition)) {
                        conditions.add(condition);
                    }
                }
            }

            if(!StringUtils.isEmpty(cellLine)) {
                IntegratedItem.setAttribute("cellLine", cellLine);
            }

            if(!StringUtils.isEmpty(IFN_gamma)) {
                IntegratedItem.setAttribute("IFN_gamma", IFN_gamma);
            }

            if(!StringUtils.isEmpty(compound)) {
                IntegratedItem.setAttribute("compound", compound);
            }

            if(!StringUtils.isEmpty(concentration) && isDouble(concentration)) {
                IntegratedItem.setAttribute("concentration", concentration);
            }

            if(!StringUtils.isEmpty(replicate)) {
                IntegratedItem.setAttribute("replicate", replicate);
            }

            if(!StringUtils.isEmpty(timepoint)) {
                IntegratedItem.setAttribute("timepoint", timepoint);
            }

            IntegratedItem.setAttribute("experiment", experimentName);

            store(IntegratedItem);
        }
    }

    private void processRNASeqExperimentGeneCount(Reader reader, String experimentName) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseTsvDelimitedReader(reader);

        String[] firstLine = (String[]) lineIter.next();
        ArrayList<String> runs = new ArrayList<String>();
        for(int i = 2; i < firstLine.length; i++) {
            String run = firstLine[i];
            runs.add(run);
        }

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String gene = line[1];

            for(int i = 2; i < line.length; i++) {
                String count = line[i];
                String runForThisItem = runs.get(i-2);

                Item IntegratedItem = createItem("RNASeqExperimentGeneCount");

                if(!gene.isEmpty()) {
                    String geneId = getGeneId(gene);

                    if (StringUtils.isEmpty(geneId)) {
                        continue;
                    }

                    IntegratedItem.setReference("gene", geneId);
                } else {
                    continue;
                }

                if(!StringUtils.isEmpty(runForThisItem)) {
                    IntegratedItem.setAttribute("run", runForThisItem);
                } else {
                    continue;
                }

                if(!StringUtils.isEmpty(count)) {
                    IntegratedItem.setAttribute("count", count);
                }

                IntegratedItem.setAttribute("experiment", experimentName);

                store(IntegratedItem);
            }
        }
    }

    private Map<String, File> readDirectoriesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            if(file.isDirectory()) {
                files.put(file.getName(), file);
            }
        }
        return files;
    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private String getGeneId(String primaryIdentifier) throws ObjectStoreException {
        try {
            String resolvedIdentifier = resolveGene(primaryIdentifier);
            if (StringUtils.isEmpty(resolvedIdentifier)) {
                return null;
            }
            String geneId = genes.get(resolvedIdentifier);
            if (geneId == null) {
                Item gene = createItem("Gene");
                gene.setAttribute("primaryIdentifier", resolvedIdentifier);
                store(gene);
                geneId = gene.getIdentifier();
                genes.put(resolvedIdentifier, geneId);
            }
            return geneId;
        } catch (Exception e) {
            return "";
        }
    }

    private String resolveGene(String identifier) {
        String id = identifier;

        if (rslv != null && rslv.hasTaxon(TAXON_ID)) {
            int resCount = rslv.countResolutions(TAXON_ID, identifier);
            if (resCount != 1) {
                LOG.info("RESOLVER: failed to resolve gene to one identifier, ignoring gene: "
                        + identifier + " count: " + resCount + " Human identifier: "
                        + rslv.resolveId(TAXON_ID, identifier));
                return null;
            }
            id = rslv.resolveId(TAXON_ID, identifier).iterator().next();
        }
        return id;
    }

    private boolean isDouble(String str) {
        try {
            double x = Double.parseDouble(str);
            if (x == (int) x)
                return false;
            return true;
        }
        catch(NumberFormatException e) {
            return false;
        }
    }
}
