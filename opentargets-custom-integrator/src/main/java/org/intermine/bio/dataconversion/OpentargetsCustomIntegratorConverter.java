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
public class OpentargetsCustomIntegratorConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "DisGeNET";
    private static final String DATA_SOURCE_NAME = "Curated gene-disease associations";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String ASSOCIATIONS_SCORE_FILE = "opentargets_evidence_storm_targets.csv";
    private static final String ASSOCIATIONS_PAPERS_FILE = "opentargets_papers_storm_targets.tsv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, Item> diseases = new HashMap<String, Item>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(OpentargetsCustomIntegratorConverter.class);

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.


    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public OpentargetsCustomIntegratorConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(TAXON_ID);
        }
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {

        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processAssociationsScores(new FileReader(files.get(ASSOCIATIONS_SCORE_FILE)));
        processAssociationsPapers(new FileReader(files.get(ASSOCIATIONS_PAPERS_FILE)));

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processAssociationsScores(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // Skip header
        lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String geneSymbol = line[0];
            String diseaseName = line[1];
            String score = line[2];
            String nrReferences = line[3];

            //Item disease = getDisease(diseaseId);
            Item disease = null;

            if (diseases.get(diseaseName) == null) {
                disease = createItem("Disease");
                disease.setAttribute("name", diseaseName);
                try {
                    store(disease);
                } catch (ObjectStoreException e) {
                    throw new RuntimeException("failed to store disease with primary identifier: " + diseaseName, e);
                }
                diseases.put(diseaseName, disease);

            } else {
                disease = getDisease(diseaseName);
            }

            String geneId = getGeneId(geneSymbol);

            if (StringUtils.isEmpty(geneId)) {
                continue;
            }

            Item interactionItem;
            interactionItem = createItem("OpenTargetsDiseaseAssociationsScore");
            interactionItem.setReference("gene", geneId);
            interactionItem.setReference("disease", disease);
            interactionItem.setAttribute("score", score);
            interactionItem.setAttribute("nrReferences", nrReferences);
            try {
                store(interactionItem);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store diseaseassociation with disease identifier: " + diseaseName, e);
            }
        }
    }

    private void processAssociationsPapers(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // Skip header
        lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String geneSymbol = line[0];
            String diseaseName = line[1];
            String pmid = line[2];
            String url = line[3];
            String doi = line[4];
            String authorString = line[5];
            String journalTitle = line[6];
            String pubYear = line[7];
            String title = line[8];
            String citedByCount = line[9];

            //Item disease = getDisease(diseaseId);
            Item disease = null;

            if (diseases.get(diseaseName) == null) {
                disease = createItem("Disease");
                disease.setAttribute("name", diseaseName);
                try {
                    store(disease);
                } catch (ObjectStoreException e) {
                    throw new RuntimeException("failed to store disease with primary identifier: " + diseaseName, e);
                }
                diseases.put(diseaseName, disease);

            } else {
                disease = getDisease(diseaseName);
            }

            String geneId = getGeneId(geneSymbol);

            if (StringUtils.isEmpty(geneId)) {
                continue;
            }

            Item interactionItem;
            interactionItem = createItem("OpenTargetsDiseaseAssociationsPaper");
            interactionItem.setReference("gene", geneId);
            interactionItem.setReference("disease", disease);
            interactionItem.setAttribute("pmid", pmid);
            interactionItem.setAttribute("url", url);
            interactionItem.setAttribute("doi", doi);
            interactionItem.setAttribute("authorString", authorString);
            interactionItem.setAttribute("journalTitle", journalTitle);
            interactionItem.setAttribute("pubYear", pubYear);
            interactionItem.setAttribute("title", title);
            interactionItem.setAttribute("citedByCount", citedByCount);
            try {
                store(interactionItem);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store diseaseassociation with disease identifier: " + diseaseName, e);
            }
        }
    }

    private String getGeneId(String primaryIdentifier) throws ObjectStoreException {
        String resolvedIdentifier = resolveGene(primaryIdentifier);
        if (StringUtils.isEmpty(resolvedIdentifier)) {
            return null;
        }
        String geneId = genes.get(resolvedIdentifier);
        if (geneId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("primaryIdentifier", resolvedIdentifier);
            gene.setReference("organism", getOrganism(TAXON_ID));
            store(gene);
            geneId = gene.getIdentifier();
            genes.put(resolvedIdentifier, geneId);
        }
        return geneId;
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

    private Item getDisease(String diseaseId) {
        Item disease = diseases.get(diseaseId);
        return disease;
    }
}
