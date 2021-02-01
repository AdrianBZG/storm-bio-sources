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
import java.util.ArrayList;
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
    private static final String DATASET_TITLE = "OpenTargets";
    private static final String DATA_SOURCE_NAME = "Gene-disease associations and their scores from OpenTargets";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String ASSOCIATIONS_SCORE_FILE = "opentargets_evidence_storm_targets.csv";
    private static final String ASSOCIATIONS_PAPERS_FILE = "opentargets_papers_storm_targets.csv";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, Item> diseases = new HashMap<String, Item>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(OpentargetsCustomIntegratorConverter.class);

    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    // Methods to integrate the data only for a list of genes
    private static final String GENE_LIST_FILE = "/data/storm/targets/storm_targets_symbols.csv";
    private ArrayList<String> processGeneList(String geneListFile) throws Exception {
        File geneListF = new File(geneListFile);

        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(new FileReader(geneListF));
        ArrayList<String> geneListArray = new ArrayList<String>();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();
            String gene = line[0];
            if(StringUtils.isEmpty(gene)) {
                continue;
            }

            String resolvedGeneIdentifier = getGeneIdentifier(gene);
            if(resolvedGeneIdentifier != null) {
                geneListArray.add(resolvedGeneIdentifier);
            }
        }

        return geneListArray;
    }

    private String getGeneIdentifier(String geneSymbol) throws ObjectStoreException {
        String resolvedIdentifier = resolveGene(geneSymbol);
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
    }
    //

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

        ArrayList<String> geneListArray = new ArrayList<String>();
        if(!StringUtils.isEmpty(GENE_LIST_FILE)) {
            geneListArray = processGeneList(GENE_LIST_FILE);
        }

        processAssociationsScores(new FileReader(files.get(ASSOCIATIONS_SCORE_FILE)), geneListArray);
        processAssociationsPapers(new FileReader(files.get(ASSOCIATIONS_PAPERS_FILE)), geneListArray);

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    private void processAssociationsScores(Reader reader, ArrayList<String> geneList) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // Skip header
        lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String geneSymbol = line[1];

            if(!geneList.isEmpty()) {
                String resolvedGene = getGeneIdentifier(geneSymbol);
                if(!geneList.contains(resolvedGene)) {
                    continue;
                }
            }

            String diseaseName = line[2];
            String assocScore = line[3];
            //String score = String.valueOf(scoreNumber);
            String rnaScore = line[4];
            String geneticScore = line[5];
            String somaticScore = line[6];
            String knownDrugScore = line[7];
            String animalModelScore = line[8];
            String affectedPathwayScore = line[9];
            String litScore = line[10];
            //String litScore = String.valueOf(litScoreNumber);
            int referencesNumber = new Double(line[11]).intValue();
            //int referencesNumberInt = (int)referencesNumber;
            String nrReferences = String.valueOf(referencesNumber);

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
            interactionItem = createItem("OpenTargetsScores");
            interactionItem.setReference("gene", geneId);
            interactionItem.setReference("disease", disease);
            interactionItem.setAttribute("score", assocScore);
            interactionItem.setAttribute("rnaScore", rnaScore);
            interactionItem.setAttribute("geneticScore", geneticScore);
            interactionItem.setAttribute("somaticScore", somaticScore);
            interactionItem.setAttribute("knownDrugScore", knownDrugScore);
            interactionItem.setAttribute("animalModelScore", animalModelScore);
            interactionItem.setAttribute("affectedPathwayScore", affectedPathwayScore);
            interactionItem.setAttribute("litScore", litScore);
            interactionItem.setAttribute("nrReferences", String.valueOf(nrReferences));
            try {
                store(interactionItem);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store diseaseassociation with disease identifier: " + diseaseName, e);
            }
        }
    }

    private void processAssociationsPapers(Reader reader, ArrayList<String> geneList) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // Skip header
        lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String geneSymbol = line[1];

            if(!geneList.isEmpty()) {
                String resolvedGene = getGeneIdentifier(geneSymbol);
                if(!geneList.contains(resolvedGene)) {
                    continue;
                }
            }

            String diseaseName = line[2];
            String pmid = line[3];
            String url = line[5];
            String doi = line[4];
            String authorString = line[6];
            String journalTitle = line[7];
            String pubYear = line[8];
            String title = line[9];
            Integer citedByCountNumber = Integer.parseInt(line[10]);
            String citedByCount = String.valueOf(citedByCountNumber);

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

            if (StringUtils.isEmpty(pmid) || StringUtils.isEmpty(url) || StringUtils.isEmpty(doi) || StringUtils.isEmpty(authorString) || StringUtils.isEmpty(journalTitle) || StringUtils.isEmpty(pubYear) || StringUtils.isEmpty(title) || StringUtils.isEmpty(citedByCount)) {
                continue;
            }

            Item interactionItem;
            interactionItem = createItem("OpenTargetsEvidences");
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
            //gene.setAttribute("symbol", primaryIdentifier);
            //gene.setReference("organism", getOrganism(TAXON_ID));
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
