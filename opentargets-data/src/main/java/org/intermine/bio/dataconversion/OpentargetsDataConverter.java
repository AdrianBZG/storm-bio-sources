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
import java.util.*;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;
import org.json.*;
import java.io.BufferedReader;
import org.json.JSONObject;
// Parsing JSONs: studytrails.com/java/json/java-org-json/


/**
 * 
 * @author
 */
public class OpentargetsDataConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "OpenTargets Associations";
    private static final String DATA_SOURCE_NAME = "Data on OpenTargets gene-disease associations";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private static final String DISEASES_FILE = "19.11_disease_list.csv";
    private static final String ASSOCIATIONS_FILE = "19.11_association_data.json";

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, Item> diseases = new HashMap<String, Item>();


    private String organismIdentifier; // Not the taxon ID. It references the object that is created into the database.

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public OpentargetsDataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * 
     *
     * {@inheritDoc}
     */
    public void process(File dataDir) throws Exception {
        Map<String, File> files = readFilesInDir(dataDir);

        organismIdentifier = getOrganism(TAXON_ID);

        processDiseases(new FileReader(files.get(DISEASES_FILE)));
        processAssociations(new FileReader(files.get(ASSOCIATIONS_FILE)));

    }

    private Map<String, File> readFilesInDir(File dir) {
        Map<String, File> files = new HashMap<String, File>();
        for (File file : dir.listFiles()) {
            files.put(file.getName(), file);
        }
        return files;
    }

    // Load the diseases in OpenTargets
    private void processDiseases(Reader reader) throws ObjectStoreException, IOException {
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(reader);
        // Skip header
        //lineIter.next();

        while (lineIter.hasNext()) {
            String[] line = (String[]) lineIter.next();

            String diseaseId = line[0];
            String diseaseName = line[1];

            Item disease = getDisease(diseaseId);

            if (disease == null) {
                disease = createItem("Disease");
                disease.setAttribute("diseaseId", diseaseId);
                disease.setAttribute("diseaseName", diseaseName);
                disease.setAttribute("diseaseType", "NA");
                disease.setAttribute("primaryIdentifier", diseaseId);
                store(disease);
            }
        }
    }

    // Load the associations in OpenTargets
    private void processAssociations(Reader reader) throws ObjectStoreException, IOException {
        try {
            try (BufferedReader br = new BufferedReader(reader)) {
                String line;
                while ((line = br.readLine()) != null) {
                    try {
                        JSONObject jsonObject = new JSONObject(line);

                        // Get the gene symbol
                        JSONObject geneSymbolJson = (JSONObject)jsonObject.get("target");
                        geneSymbolJson = (JSONObject)geneSymbolJson.get("gene_info");
                        String geneSymbol = (String)geneSymbolJson.get("symbol");

                        // Get the overall association score
                        JSONObject overallAssociationScoreJson = (JSONObject)jsonObject.get("association_score");
                        String overallAssociationScore = (String)overallAssociationScoreJson.get("overall");

                        // Scores for each data type
                        JSONObject overallAssociationScoreDatatypesJson = (JSONObject)overallAssociationScoreJson.get("datatypes");
                        String literatureScore = (String)overallAssociationScoreDatatypesJson.get("literature");
                        String rnaExpressionScore = (String)overallAssociationScoreDatatypesJson.get("rna_expression");
                        String geneticAssociationScore = (String)overallAssociationScoreDatatypesJson.get("genetic_association");
                        String somaticMutationScore = (String)overallAssociationScoreDatatypesJson.get("somatic_mutation");
                        String knownDrugScore = (String)overallAssociationScoreDatatypesJson.get("known_drug");
                        String animalModelScore = (String)overallAssociationScoreDatatypesJson.get("animal_model");
                        String affectedPathwayScore = (String)overallAssociationScoreDatatypesJson.get("affected_pathway");

                        // Counts for each data type
                        JSONObject evidenceCountJson = (JSONObject)jsonObject.get("evidence_count");
                        evidenceCountJson = (JSONObject)evidenceCountJson.get("datayypes");
                        String literatureCount = (String)evidenceCountJson.get("literature");
                        String rnaExpressionCount = (String)evidenceCountJson.get("rna_expression");
                        String geneticAssociationCount = (String)evidenceCountJson.get("genetic_association");
                        String somaticMutationCount = (String)evidenceCountJson.get("somatic_mutation");
                        String knownDrugCount = (String)evidenceCountJson.get("known_drug");
                        String animalModelCount = (String)evidenceCountJson.get("animal_model");
                        String affectedPathwayCount = (String)evidenceCountJson.get("affected_pathway");

                        // Get the disease ID for reference
                        JSONObject diseaseIdJson = (JSONObject)jsonObject.get("disease");
                        String diseaseId = (String)diseaseIdJson.get("id");

                        // Get the therapeutic areas and create a relation from disease id - therapeutic area for each one of them
                        // Have to check that it wasn't added yet
                        JSONObject therapeuticAreasJson = (JSONObject)diseaseIdJson.get("efo_info");
                        therapeuticAreasJson = (JSONObject)therapeuticAreasJson.get("therapeutic_area");
                        JSONArray therapeuticAreas = (JSONArray)therapeuticAreasJson.get("labels");

                        List<Object> therapeuticAreasList = toList(therapeuticAreas);
                        for(Object element : therapeuticAreasList) {
                            String therapeuticAreaString = element.toString();
                            // Check that the diseaseId in the hashmap doesn't have this therapeutic area in its value (list of areas)

                            // Else create the disease-therapeutic area relation
                        }

                    } catch (JSONException err) {
                        throw new RuntimeException("Failed to read the following JSON in OpenTargets associations: " + line, err);
                    }
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Failed to process OpenTargets associations: " + e);
        }
    }

    private Map<String, Object> toMap(JSONObject object) throws JSONException {
        Map<String, Object> map = new HashMap<String, Object>();

        Iterator<String> keysItr = object.keys();
        while(keysItr.hasNext()) {
            String key = keysItr.next();
            Object value = object.get(key);

            if(value instanceof JSONArray) {
                value = toList((JSONArray) value);
            }

            else if(value instanceof JSONObject) {
                value = toMap((JSONObject) value);
            }
            map.put(key, value);
        }
        return map;
    }


    public List<Object> toList(JSONArray array) throws JSONException {
        List<Object> list = new ArrayList<Object>();
        for(int i = 0; i < array.length(); i++) {
            Object value = array.get(i);
            if(value instanceof JSONArray) {
                value = toList((JSONArray) value);
            }

            else if(value instanceof JSONObject) {
                value = toMap((JSONObject) value);
            }
            list.add(value);
        }
        return list;
    }

    public String getGene(String identifier) {
        String refId = genes.get(identifier);
        if (refId == null) {
            Item gene = createItem("Gene");
            gene.setAttribute("symbol", identifier);
            gene.setReference("organism", getOrganism(TAXON_ID));
            try {
                store(gene);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store gene with primary identifier: " + identifier, e);
            }
            refId = gene.getIdentifier();
            genes.put(identifier, refId);
        }
        return refId;
    }

    private Item getDisease(String diseaseId) {
        Item disease = diseases.get(diseaseId);
        return disease;
    }
}
