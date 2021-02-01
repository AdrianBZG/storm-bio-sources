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
    private static final String DATASET_TITLE = "STORM Correlations Analyses";
    private static final String DATA_SOURCE_NAME = "Results for the 2020 correlations analyses on STORM Targets";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private Map<String, String> genes = new HashMap<String, String>();

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
            String folderName = entry.getKey();
            Map<String, File> files = readFilesInDir(theFolder);
            for (Map.Entry<String, File> fileEntry : files.entrySet()) {
                processCellLineRNASeqExperiment(new FileReader(fileEntry.getValue()), folderName);
                break;
            }
        }

    }

    private void processCellLineRNASeqExperiment(Reader reader, String type) throws ObjectStoreException, IOException {

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
                //gene.setAttribute("symbol", primaryIdentifier);
                //gene.setReference("organism", getOrganism(TAXON_ID));
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
}
