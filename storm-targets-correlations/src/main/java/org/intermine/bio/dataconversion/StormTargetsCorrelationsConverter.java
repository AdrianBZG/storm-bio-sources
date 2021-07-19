package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2018 FlyMine
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
public class StormTargetsCorrelationsConverter extends BioDirectoryConverter
{
    private static final String DATASET_TITLE = "STORM Correlations Analyses";
    private static final String DATA_SOURCE_NAME = "Results for the 2020 correlations analyses on STORM Targets";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormTargetsCorrelationsConverter.class);

    private String organismIdentifier;
    public StormTargetsCorrelationsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    private boolean isDouble(String str) {
        try {
            // check if it can be parsed as any double
            double x = Double.parseDouble(str);
            // check if the double can be converted without loss to an int
            if (x == (int) x)
                // if yes, this is an int, thus return false
                return false;
            // otherwise, this cannot be converted to an int (e.g. "1.2")
            return true;
            // short version: return x != (int) x;
        }
        catch(NumberFormatException e) {
            return false;
        }
    }

    public void process(File dataDir) throws Exception {
        organismIdentifier = getOrganism(TAXON_ID);

        Map<String, File> directories = readDirectoriesInDir(dataDir);

        for (Map.Entry<String, File> entry : directories.entrySet()) {
            File theFolder = entry.getValue();
            String folderName = entry.getKey();
            Map<String, File> files = readFilesInDir(theFolder);
            for (Map.Entry<String, File> fileEntry : files.entrySet()) {
                processTargetsCorrelations(new FileReader(fileEntry.getValue()), folderName);
            }
        }

    }

    private void processTargetsCorrelations(Reader reader, String type) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // header
        String[] firstLine = (String[]) lineIter.next();

        //lineIter.next();
        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String gene1 = line[0].split(" ")[0].trim();
            String gene2 = line[1].split(" ")[0].trim();
            String nsize = line[2];
            String correlation = line[3];
            String pvalue = line[4];
            String fdr = line[5];

            Item IntegratedItem;


            IntegratedItem = createItem("STORMTargetCorrelations");

            if(!gene1.isEmpty()) {
                if(unresolvableGenes.get(gene1) != null) {
                    continue;
                }
                String geneId = getGeneId(gene1);
                if(geneId == null) {
                    continue;
                }

                IntegratedItem.setReference("gene1", geneId);
            } else {
                continue;
            }


            if(!gene2.isEmpty()) {
                if(unresolvableGenes.get(gene2) != null) {
                    continue;
                }
                String geneId = getGeneId(gene2);
                if(geneId == null) {
                    continue;
                }

                IntegratedItem.setReference("gene2", geneId);
            } else {
                continue;
            }


            if(!StringUtils.isEmpty(type)) {
                IntegratedItem.setAttribute("experimentType", type);
            } else {
                continue;
            }

            if(!StringUtils.isEmpty(nsize)) {
                IntegratedItem.setAttribute("nsize", nsize);
            }

            if(!StringUtils.isEmpty(correlation) && isDouble(correlation)) {
                IntegratedItem.setAttribute("correlation", correlation);
            }

            if(!StringUtils.isEmpty(pvalue) && isDouble(pvalue)) {
                IntegratedItem.setAttribute("pvalue", pvalue);
            }

            if(!StringUtils.isEmpty(fdr) && isDouble(fdr)) {
                IntegratedItem.setAttribute("fdr", fdr);
            }


            store(IntegratedItem);
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

    private String getGeneId(String identifier) throws ObjectStoreException {
        String geneId = null;
        try {
            String resolvedIdentifier = resolveGene(identifier);
            if(resolvedIdentifier != null) {
                geneId = genes.get(resolvedIdentifier);
                if (geneId == null) {
                    Item gene = createItem("Gene");
                    gene.setAttribute("primaryIdentifier", resolvedIdentifier);
                    store(gene);
                    geneId = gene.getIdentifier();
                    genes.put(resolvedIdentifier, geneId);
                }
                return geneId;
            } else {
                return resolvedIdentifier;
            }
        } catch (Exception e) {
            LOG.info("getGeneId: failed to resolve gene: " + identifier);
            return null;
        }
    }

    private String resolveGene(String identifier) {
        String id = null;

        if(resolvedGenes.get(identifier) != null) {
            id = resolvedGenes.get(identifier);
        } else {
            if (rslv != null && rslv.hasTaxon(TAXON_ID)) {
                int resCount = rslv.countResolutions(TAXON_ID, identifier);
                if (resCount != 1) {
                    unresolvableGenes.put(identifier, identifier);
                    return null;
                }
                id = rslv.resolveId(TAXON_ID, identifier).iterator().next();
                resolvedGenes.put(identifier, id);
            }
        }
        return id;
    }
}
