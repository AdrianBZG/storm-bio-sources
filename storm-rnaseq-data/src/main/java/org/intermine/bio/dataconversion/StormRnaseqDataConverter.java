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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
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
    private String configFile = null;

    private Map<String, String> genes = new HashMap<String, String>();
    private ArrayList<String> conditions = new ArrayList<String>();

    private Map<String, String> experiments = new HashMap<String, String>();
    private Map<String, String> experimentsSampleInfoPath = new HashMap<String, String>();
    private Map<String, String> experimentsGeneCountsPath = new HashMap<String, String>();
    private Map<String, String> experimentsDESEQ2Path = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormRnaseqDataConverter.class);

    private String organismIdentifier;
    public StormRnaseqDataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }


    public void setConfigFile(String configFile) {
        this.configFile = configFile;
    }

    public void process(File dataDir) throws Exception {
        organismIdentifier = getOrganism(TAXON_ID);

        Map<String, File> directories = readDirectoriesInDir(dataDir);

        if (this.configFile == null) {
            throw new BuildException("configFile attribute is not set");
        }

        // Load the config file
        processConfigFile(new FileReader(dataDir + "/config.csv"));
        //

        for (Map.Entry<String, File> entry : directories.entrySet()) {
            File theFolder = entry.getValue();
            String experimentFolderName = entry.getKey();
            if(experiments.get(experimentFolderName) != null && experimentsSampleInfoPath.get(experimentFolderName) != null && experimentsGeneCountsPath.get(experimentFolderName) != null) {
                Map<String, File> files = readFilesInDir(theFolder);
                for (Map.Entry<String, File> fileEntry : files.entrySet()) {
                    String cellLineFolderName = fileEntry.getKey();
                    String experimentName = experimentFolderName + "-" + cellLineFolderName;

                    String rootFolder = theFolder.getAbsolutePath();

                    // Sample info
                    String sampleInfoRelativePath = experimentsSampleInfoPath.get(experimentFolderName).replace("FOLDERNAME", cellLineFolderName);
                    String sampleInfoPath = rootFolder.concat(sampleInfoRelativePath);
                    //processRNASeqExperimentSampleInfo(new FileReader(sampleInfoPath), experimentName);

                    // Gene counts
                    String geneCountsRelativePath = experimentsGeneCountsPath.get(experimentFolderName).replace("FOLDERNAME", cellLineFolderName);
                    String countsPath = rootFolder.concat(geneCountsRelativePath);
                    //processRNASeqExperimentGeneCount(new FileReader(countsPath), experimentName);

                    // DESeq2
                    String DESEQ2RelativePath = experimentsDESEQ2Path.get(experimentFolderName);
                    String DESEQ2Path = rootFolder.concat("/" + DESEQ2RelativePath);
                    processRNASeqExperimentDESEQ2(new FileReader(DESEQ2Path), experimentName);
                }
            }
        }

    }

    private void processConfigFile(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        String[] firstLine = (String[]) lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            String experimentFolder = line[0];
            String experimentSampleInfoPath = line[1];
            String experimentGeneCountsPath = line[2];
            String experimentDESEQ2Path = line[3];

            experiments.put(experimentFolder, experimentFolder);
            experimentsSampleInfoPath.put(experimentFolder, experimentSampleInfoPath);
            experimentsGeneCountsPath.put(experimentFolder, experimentGeneCountsPath);
            experimentsDESEQ2Path.put(experimentFolder, experimentDESEQ2Path);
        }
    }

    private void processRNASeqExperimentDESEQ2(Reader reader, String experimentName) throws ObjectStoreException, IOException {
        //
    }

    /*private void processRNASeqExperimentSampleInfo(Reader reader, String experimentName) throws ObjectStoreException, IOException {
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
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(reader);

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
    }*/

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
