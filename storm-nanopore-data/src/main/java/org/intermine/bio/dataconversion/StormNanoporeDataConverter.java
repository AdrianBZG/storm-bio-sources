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
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

import org.json.*;
import java.io.BufferedReader;
import org.json.JSONObject;

/**
 * 
 * @author
 */
public class StormNanoporeDataConverter extends BioDirectoryConverter
{
    //
    private static final String DATASET_TITLE = "STORM Nanopore Data";
    private static final String DATA_SOURCE_NAME = "STORM Nanopore Data";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private Map<String, String> genes = new HashMap<String, String>(); 
    private Map<String, Item> experiments = new HashMap<>();

    protected IdResolver rslv;
    private static final Logger LOG = Logger.getLogger(StormNanoporeDataConverter.class);

    private String organismIdentifier;

    public StormNanoporeDataConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    public void process(File dataDir) throws Exception {
        organismIdentifier = getOrganism(TAXON_ID);

        Map<String, File> directories = readDirectoriesInDir(dataDir);

        // Get all JSON config files in the directory and process one by one
        File[] configFilesArray = dataDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".json");
            }
        });

        ArrayList<File> configFiles = new ArrayList<File>();
        for(int i = 0; i < configFilesArray.length; i++) {
            configFiles.add(configFilesArray[i]);
        }

        for (File configFile : configFiles) {
            processConfigFile(new FileReader(configFile), dataDir);
        }
    }

    private void processConfigFile(Reader reader, File dataDir) throws ObjectStoreException, IOException {
        try (BufferedReader br = new BufferedReader(reader)) {
            String line;
            while ((line = br.readLine()) != null) {
                try {
                    JSONObject jsonObject = new JSONObject(line);

                    // Get the experiment key
                    JSONObject experimentJson = jsonObject.getJSONObject("experiment");

                    // Get the experiment metadata
                    String experimentName = experimentJson.getString("name");
                    String experimentShortName = experimentJson.getString("short name");
                    String experimentProject = experimentJson.getString("project");
                    String experimentContactPerson = experimentJson.getString("contact person");
                    String experimentDate = experimentJson.getString("date");
                    String experimentSequencing = experimentJson.getString("sequencing");
                    String experimentProvider = experimentJson.getString("provider");
                    String experimentDotmaticsReference = experimentJson.getString("Dotmatics reference");

                    // Save the item
                    Item ExperimentMetadataItem = createItem("NanoporeExperimentMetadata");                    
                    ExperimentMetadataItem.setAttribute("name", experimentName);
                    ExperimentMetadataItem.setAttribute("shortName", experimentShortName);
                    ExperimentMetadataItem.setAttribute("project", experimentProject);
                    ExperimentMetadataItem.setAttribute("contactPerson", experimentContactPerson);
                    ExperimentMetadataItem.setAttribute("date", experimentDate);
                    ExperimentMetadataItem.setAttribute("sequencing", experimentSequencing);
                    ExperimentMetadataItem.setAttribute("provider", experimentProvider);
                    ExperimentMetadataItem.setAttribute("dotmaticsReference", experimentDotmaticsReference);

                    String experimentKey = experimentShortName;
                    if (!experiments.containsKey(experimentKey)) {
                        experiments.put(experimentKey, ExperimentMetadataItem);
                    }

                    // Process materials
                    JSONObject materialsJson = jsonObject.getJSONObject("materials");
                    processExperimentMaterials(materialsJson, experimentShortName);

                    // Process treatments
                    JSONObject treatmentsJson = jsonObject.getJSONObject("treatments");
                    processExperimentTreatments(treatmentsJson, experimentShortName);

                    // Process each comparison individually                    
                    JSONArray experimentComparisons = (JSONArray)jsonObject.get("comparisons");
                    for(int i = 0; i < experimentComparisons.length(); i++) {
                        JSONObject comparison = experimentComparisons.getJSONObject(i);
                        JSONObject treatmentObject = comparison.getJSONObject("treatment");
                        JSONObject controlObject = comparison.getJSONObject("control");
                        String treatmentName = treatmentObject.getString("name");
                        String controlName = controlObject.getString("name");
                        
                        String comparisonName = treatmentName + "_vs_" + controlName;
                        Map<String, File> filesInDir = readFilesInDir(new File(dataDir.getAbsolutePath() + "/" + experimentShortName + "/" + comparisonName));

                        // Nanocompore output
                        String nanocomporeResultsFile = "out_nanocompore_results.tsv";
                        if(filesInDir.get(nanocomporeResultsFile) != null) {
                            File NanocomporeResultsFile = filesInDir.get(nanocomporeResultsFile);
                            processNanocomporeResults(NanocomporeResultsFile, experimentShortName);
                        } else {
                            LOG.info("Failed to find NanocomporeResultsFile file: " + nanocomporeResultsFile);
                            continue;
                        }

                        // Nanocompore Insig Results
                        String nanocomporeInsigFile = "insigResults.csv";
                        if(filesInDir.get(nanocomporeInsigFile) != null) {
                            File NanocomporeInsigFile = filesInDir.get(nanocomporeInsigFile);
                            processNanocomporeInsigFile(NanocomporeInsigFile, experimentShortName);
                        } else {
                            LOG.info("Failed to find nanocomporeInsigFile file: " + nanocomporeInsigFile);
                            continue;
                        }

                        // Nanocompore Sig Results
                        String nanocomporeSigFile = "sigResultsOrderedByLFC.csv";
                        if(filesInDir.get(nanocomporeSigFile) != null) {
                            File NanocomporeSigFile = filesInDir.get(nanocomporeSigFile);
                            processNanocomporeSigFile(NanocomporeSigFile, experimentShortName);
                        } else {
                            LOG.info("Failed to find nanocomporeSigFile file: " + nanocomporeSigFile);
                            continue;
                        }

                        // NanoporeExperimentTranscriptCounts
                        String nanoporeExperimentTranscriptCountsFile = "masterTranscriptCounts.txt";
                        if(filesInDir.get(nanoporeExperimentTranscriptCountsFile) != null) {
                            File NanoporeExperimentTranscriptCountsFile = filesInDir.get(nanoporeExperimentTranscriptCountsFile);
                            processExperimentTranscriptCount(NanoporeExperimentTranscriptCountsFile, experimentShortName);
                        } else {
                            LOG.info("Failed to find nanoporeExperimentTranscriptCountsFile file: " + nanoporeExperimentTranscriptCountsFile);
                            continue;
                        }
                    }
                    
                    storeExperiment(ExperimentMetadataItem, experimentShortName);

                } catch (JSONException err) {
                    throw new RuntimeException("Failed to read the following JSON: " + line, err);
                }
            }
        }
    }

    private void storeExperiment(Item ExperimentMetadataItem, String experimentKey) {
        try {
            store(ExperimentMetadataItem);
        } catch (Exception e) {
            throw new RuntimeException("Error storing StormNanoporeExperiment ", e);
        }
    }

    private void processExperimentTranscriptCount(File NanoporeExperimentTranscriptCountsFile, String experimentShortName) throws ObjectStoreException, IOException {
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        String fileAbsPath = NanoporeExperimentTranscriptCountsFile.getAbsolutePath();
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] firstLine = (String[]) lineIter.next();

        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            try {
                String transcript = line[0];
                String BC1 = line[1];
                String BC2 = line[2];
                String BC3 = line[3];
                String BC4 = line[4];

                Item IntegratedItem = createItem("NanoporeExperimentTranscriptCounts");
                if(!transcript.isEmpty()) {
                    IntegratedItem.setAttribute("transcript", transcript);
                } else {
                    continue;
                }

                if(!StringUtils.isEmpty(BC1) && isDouble(BC1)) {
                    IntegratedItem.setAttribute("BC1", BC1);
                }

                if(!StringUtils.isEmpty(BC2) && isDouble(BC2)) {
                    IntegratedItem.setAttribute("BC2", BC2);
                }

                if(!StringUtils.isEmpty(BC3) && isDouble(BC3)) {
                    IntegratedItem.setAttribute("BC3", BC3);
                }

                if(!StringUtils.isEmpty(BC4) && isDouble(BC4)) {
                    IntegratedItem.setAttribute("BC4", BC4);
                }

                IntegratedItem.setReference("experiment", ExperimentMetadataItem);

                store(IntegratedItem);
                
            } catch (Exception e) {
                LOG.info("Exception in processExperimentTranscriptCount with: " + line + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processNanocomporeResults(File NanocomporeResultsFile, String experimentShortName) throws ObjectStoreException, IOException {
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        String fileAbsPath = NanocomporeResultsFile.getAbsolutePath();
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] firstLine = (String[]) lineIter.next();

        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            try {
                /*
                <reference name="experiment" referenced-type="NanoporeExperimentMetadata" reverse-reference="nanoPoreExperimentNanocompores"/>
                <attribute name="pos" type="java.lang.Double"/>
                <attribute name="ref_id" type="java.lang.String"/>
                <attribute name="ref_kmer" type="java.lang.String"/>
                <attribute name="GMM_anova_pvalue" type="java.lang.Double"/>
                <attribute name="GMM_logit_pvalue" type="java.lang.Double"/>
                <attribute name="KS_dwell_pvalue" type="java.lang.Double"/>
                <attribute name="KS_intensity_pvalue" type="java.lang.Double"/>
                <attribute name="GMM_cov_type" type="java.lang.String"/>
                <attribute name="GMM_n_clust" type="java.lang.Double"/>
                <attribute name="cluster_counts" type="java.lang.String"/>
                <attribute name="Anova_delta_logit" type="java.lang.Double"/>
                <attribute name="Logit_LOR" type="java.lang.Double"/>  
                */
                String pos = line[0];
                String ref_id = line[3];
                String ref_kmer = line[5];
                String GMM_anova_pvalue = line[6];
                String GMM_logit_pvalue = line[7];
                String KS_dwell_pvalue = line[8];
                String KS_intensity_pvalue = line[9];
                String GMM_cov_type = line[10];
                String GMM_n_clust = line[11];
                String cluster_counts = line[12];
                String Anova_delta_logit = line[13];
                String Logit_LOR = line[14];

                Item IntegratedItem = createItem("NanoporeExperimentNanocompore");

                if(!StringUtils.isEmpty(pos) && isDouble(pos)) {
                    IntegratedItem.setAttribute("pos", pos);
                }

                if(!ref_id.isEmpty()) {
                    IntegratedItem.setAttribute("ref_id", ref_id);
                }

                if(!ref_kmer.isEmpty()) {
                    IntegratedItem.setAttribute("ref_kmer", ref_kmer);
                }

                if(!GMM_cov_type.isEmpty()) {
                    IntegratedItem.setAttribute("GMM_cov_type", GMM_cov_type);
                }

                if(!cluster_counts.isEmpty()) {
                    IntegratedItem.setAttribute("cluster_counts", cluster_counts);
                }

                if(!StringUtils.isEmpty(GMM_anova_pvalue) && isDouble(GMM_anova_pvalue)) {
                    IntegratedItem.setAttribute("GMM_anova_pvalue", GMM_anova_pvalue);
                }

                if(!StringUtils.isEmpty(GMM_logit_pvalue) && isDouble(GMM_logit_pvalue)) {
                    IntegratedItem.setAttribute("GMM_logit_pvalue", GMM_logit_pvalue);
                }

                if(!StringUtils.isEmpty(KS_dwell_pvalue) && isDouble(KS_dwell_pvalue)) {
                    IntegratedItem.setAttribute("KS_dwell_pvalue", KS_dwell_pvalue);
                }

                if(!StringUtils.isEmpty(KS_intensity_pvalue) && isDouble(KS_intensity_pvalue)) {
                    IntegratedItem.setAttribute("KS_intensity_pvalue", KS_intensity_pvalue);
                }

                if(!StringUtils.isEmpty(GMM_n_clust) && isDouble(GMM_n_clust)) {
                    IntegratedItem.setAttribute("GMM_n_clust", GMM_n_clust);
                }

                if(!StringUtils.isEmpty(Anova_delta_logit) && isDouble(Anova_delta_logit)) {
                    IntegratedItem.setAttribute("Anova_delta_logit", Anova_delta_logit);
                }

                if(!StringUtils.isEmpty(Logit_LOR) && isDouble(Logit_LOR)) {
                    IntegratedItem.setAttribute("Logit_LOR", Logit_LOR);
                }
                
                IntegratedItem.setReference("experiment", ExperimentMetadataItem);               

                store(IntegratedItem);
                
            } catch (Exception e) {
                LOG.info("Exception in processNanocomporeResults with: " + line + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processNanocomporeInsigFile(File NanocomporeInsigFile, String experimentShortName) throws ObjectStoreException, IOException {
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        String fileAbsPath = NanocomporeInsigFile.getAbsolutePath();
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(new FileReader(fileAbsPath));
        String[] firstLine = (String[]) lineIter.next();

        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            try {
                String transcript = line[1];
                String baseMean = line[2];
                String log2FoldChange = line[3];
                String lfcSE = line[4];
                String stat = line[5];
                String pvalue = line[6];
                String padj = line[7];
                String FoldChange = line[8];
                String BC1 = line[9];
                String BC2 = line[10];
                String BC3 = line[11];
                String BC4 = line[12];

                Item IntegratedItem = createItem("NanoporeExperimentInsigResults");
                if(!transcript.isEmpty()) {
                    IntegratedItem.setAttribute("transcript", transcript);
                } else {
                    continue;
                }

                if(!StringUtils.isEmpty(baseMean) && isDouble(baseMean)) {
                    IntegratedItem.setAttribute("baseMean", baseMean);
                }

                if(!StringUtils.isEmpty(log2FoldChange) && isDouble(log2FoldChange)) {
                    IntegratedItem.setAttribute("log2FoldChange", log2FoldChange);
                }

                if(!StringUtils.isEmpty(lfcSE) && isDouble(lfcSE)) {
                    IntegratedItem.setAttribute("lfcSE", lfcSE);
                }

                if(!StringUtils.isEmpty(stat) && isDouble(stat)) {
                    IntegratedItem.setAttribute("stat", stat);
                }

                if(!StringUtils.isEmpty(pvalue) && isDouble(pvalue)) {
                    IntegratedItem.setAttribute("pvalue", pvalue);
                }

                if(!StringUtils.isEmpty(padj) && isDouble(padj)) {
                    IntegratedItem.setAttribute("padj", padj);
                }

                if(!StringUtils.isEmpty(FoldChange) && isDouble(FoldChange)) {
                    IntegratedItem.setAttribute("FoldChange", FoldChange);
                }

                if(!StringUtils.isEmpty(BC1) && isDouble(BC1)) {
                    IntegratedItem.setAttribute("BC1", BC1);
                }

                if(!StringUtils.isEmpty(BC2) && isDouble(BC2)) {
                    IntegratedItem.setAttribute("BC2", BC2);
                }

                if(!StringUtils.isEmpty(BC3) && isDouble(BC3)) {
                    IntegratedItem.setAttribute("BC3", BC3);
                }

                if(!StringUtils.isEmpty(BC4) && isDouble(BC4)) {
                    IntegratedItem.setAttribute("BC4", BC4);
                }                

                IntegratedItem.setReference("experiment", ExperimentMetadataItem);

                store(IntegratedItem);
                
            } catch (Exception e) {
                LOG.info("Exception in processNanocomporeInsigFile with: " + line + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processNanocomporeSigFile(File NanocomporeSigFile, String experimentShortName) throws ObjectStoreException, IOException {
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        String fileAbsPath = NanocomporeSigFile.getAbsolutePath();
        Iterator<?> lineIter = FormattedTextParser.parseCsvDelimitedReader(new FileReader(fileAbsPath));
        String[] firstLine = (String[]) lineIter.next();

        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            try {
                String transcript = line[1];
                String baseMean = line[2];
                String log2FoldChange = line[3];
                String lfcSE = line[4];
                String stat = line[5];
                String pvalue = line[6];
                String padj = line[7];
                String FoldChange = line[8];
                String BC1 = line[9];
                String BC2 = line[10];
                String BC3 = line[11];
                String BC4 = line[12];

                Item IntegratedItem = createItem("NanoporeExperimentSigResults");
                if(!transcript.isEmpty()) {
                    IntegratedItem.setAttribute("transcript", transcript);
                } else {
                    continue;
                }

                if(!StringUtils.isEmpty(baseMean) && isDouble(baseMean)) {
                    IntegratedItem.setAttribute("baseMean", baseMean);
                }

                if(!StringUtils.isEmpty(log2FoldChange) && isDouble(log2FoldChange)) {
                    IntegratedItem.setAttribute("log2FoldChange", log2FoldChange);
                }

                if(!StringUtils.isEmpty(lfcSE) && isDouble(lfcSE)) {
                    IntegratedItem.setAttribute("lfcSE", lfcSE);
                }

                if(!StringUtils.isEmpty(stat) && isDouble(stat)) {
                    IntegratedItem.setAttribute("stat", stat);
                }

                if(!StringUtils.isEmpty(pvalue) && isDouble(pvalue)) {
                    IntegratedItem.setAttribute("pvalue", pvalue);
                }

                if(!StringUtils.isEmpty(padj) && isDouble(padj)) {
                    IntegratedItem.setAttribute("padj", padj);
                }

                if(!StringUtils.isEmpty(FoldChange) && isDouble(FoldChange)) {
                    IntegratedItem.setAttribute("FoldChange", FoldChange);
                }

                if(!StringUtils.isEmpty(BC1) && isDouble(BC1)) {
                    IntegratedItem.setAttribute("BC1", BC1);
                }

                if(!StringUtils.isEmpty(BC2) && isDouble(BC2)) {
                    IntegratedItem.setAttribute("BC2", BC2);
                }

                if(!StringUtils.isEmpty(BC3) && isDouble(BC3)) {
                    IntegratedItem.setAttribute("BC3", BC3);
                }

                if(!StringUtils.isEmpty(BC4) && isDouble(BC4)) {
                    IntegratedItem.setAttribute("BC4", BC4);
                }                

                IntegratedItem.setReference("experiment", ExperimentMetadataItem);

                store(IntegratedItem);
                
            } catch (Exception e) {
                LOG.info("Exception in processNanocomporeSigFile with: " + line + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processExperimentMaterials(JSONObject materialsJson, String experimentShortName) {
        // Iterate over each material
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        Iterator<String> keys = materialsJson.keys();

        while(keys.hasNext()) {
            String key = keys.next();
            try {
                // Now we have the material
                String materialName = key;

                // There should only be one key under this
                JSONObject materialTypeKeysJSON = materialsJson.getJSONObject(key);
                ArrayList<String> materialTypeKeys = new ArrayList<String>();

                Iterator<?> iterator = materialTypeKeysJSON.keys();
                while (iterator.hasNext()) {
                    Object keyObj = iterator.next();
                    materialTypeKeys.add(key.toString());
                }

                if(materialTypeKeys.size() != 1) {
                    LOG.info("Material did not had one key under it: " + materialName);
                }

                String materialType = materialTypeKeys.get(0);
                
                // Get the material object
                JSONObject materialObject = materialTypeKeysJSON.getJSONObject(materialType);

                switch(materialType) {
                    case "cell line":
                    {
                        String cellLineName = materialObject.getString("name");
                        String cellLineTissue = materialObject.getString("tissue");
                        String cellLineSpecies = materialObject.getString("species");

                        // Save the item
                        Item MaterialMetadataItem = createItem("NanoporeExperimentMaterial");

                        MaterialMetadataItem.setAttribute("materialType", materialType);
                        MaterialMetadataItem.setAttribute("name", cellLineName);
                        MaterialMetadataItem.setAttribute("tissue", cellLineTissue);
                        MaterialMetadataItem.setAttribute("species", cellLineSpecies);

                        MaterialMetadataItem.setReference("experiment", ExperimentMetadataItem);

                        store(MaterialMetadataItem);
                    }
                    case "tumour":
                    {
                        String tumourPrimaryDisease = materialObject.getString("primary disease");
                        String tumourDiseaseSubtype = materialObject.getString("disease subtype");
                        String tumourTissue = materialObject.getString("tissue");
                        String tumourSpecies = materialObject.getString("species");

                        // Save the item
                        Item MaterialMetadataItem = createItem("NanoporeExperimentMaterial");

                        MaterialMetadataItem.setAttribute("materialType", materialType);           
                        MaterialMetadataItem.setAttribute("primaryDisease", tumourPrimaryDisease);
                        MaterialMetadataItem.setAttribute("diseaseSubtype", tumourDiseaseSubtype);
                        MaterialMetadataItem.setAttribute("tissue", tumourTissue);
                        MaterialMetadataItem.setAttribute("species", tumourSpecies);

                        MaterialMetadataItem.setReference("experiment", ExperimentMetadataItem);

                        store(MaterialMetadataItem);
                    }
                    case "tissue":
                    {
                        String tissueTissue = materialObject.getString("tissue");
                        String tissueSpecies = materialObject.getString("species");

                        // Save the item
                        Item MaterialMetadataItem = createItem("NanoporeExperimentMaterial");

                        MaterialMetadataItem.setAttribute("materialType", materialType);
                        MaterialMetadataItem.setAttribute("tissue", tissueTissue);
                        MaterialMetadataItem.setAttribute("species", tissueSpecies);

                        MaterialMetadataItem.setReference("experiment", ExperimentMetadataItem);

                        store(MaterialMetadataItem);
                    }
                }
            } catch (Exception e) {
                LOG.info("Exception in processExperimentMaterials with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processExperimentTreatments(JSONObject treatmentsJson, String experimentShortName) {
        // Iterate over each material
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        Iterator<String> keys = treatmentsJson.keys();

        while(keys.hasNext()) {
            String key = keys.next();
            try {
                // Now we have the material
                String treatmentName = key;

                // There should only be one key under this
                JSONObject treatmentTypeKeysJSON = treatmentsJson.getJSONObject(key);
                ArrayList<String> treatmentTypeKeys = new ArrayList<String>();

                Iterator<?> iterator = treatmentTypeKeysJSON.keys();
                while (iterator.hasNext()) {
                    Object keyObj = iterator.next();
                    treatmentTypeKeys.add(key.toString());
                }

                if(treatmentTypeKeys.size() != 1) {
                    LOG.info("Treatment did not had one key under it: " + treatmentName);
                }

                String treatmentType = treatmentTypeKeys.get(0);
                
                // Get the material object
                JSONObject treatmentObject = treatmentTypeKeysJSON.getJSONObject(treatmentType);

                switch(treatmentType) {
                    case "inhibitor":
                    {
                        String inhibitorName = treatmentObject.getString("name");
                        String targetGene = treatmentObject.getString("target gene");
                        String dotmaticsReference = treatmentObject.getString("Dotmatics reference");
                        String dose = treatmentObject.getString("dose");
                        String timePoint = treatmentObject.getString("time point");                            

                        // Save the item
                        Item TreatmentMetadataItem = createItem("NanoporeExperimentTreatment");

                        TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);
                        TreatmentMetadataItem.setAttribute("name", key);
                        TreatmentMetadataItem.setAttribute("targetGene", targetGene);
                        TreatmentMetadataItem.setAttribute("dotmaticsReference", dotmaticsReference);

                        if(!StringUtils.isEmpty(dose) && isDouble(dose)) {
                            TreatmentMetadataItem.setAttribute("dose_concentration", dose);
                        }

                        TreatmentMetadataItem.setAttribute("timePoint", timePoint);

                        TreatmentMetadataItem.setReference("experiment", ExperimentMetadataItem);

                        store(TreatmentMetadataItem);
                    }
                    case "knock-down":
                    {
                        String inhibitorName = treatmentObject.getString("name");
                        String targetGene = treatmentObject.getString("target gene");
                        String concentration = treatmentObject.getString("concentration");
                        String type = treatmentObject.getString("type");
                        String timePoint = treatmentObject.getString("time point");                            

                        // Save the item
                        Item TreatmentMetadataItem = createItem("NanoporeExperimentTreatment");

                        TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);
                        TreatmentMetadataItem.setAttribute("name", key);
                        TreatmentMetadataItem.setAttribute("targetGene", targetGene);

                        if(!StringUtils.isEmpty(concentration) && isDouble(concentration)) {
                            TreatmentMetadataItem.setAttribute("dose_concentration", concentration);
                        }

                        TreatmentMetadataItem.setAttribute("type", type);
                        TreatmentMetadataItem.setAttribute("timePoint", timePoint);

                        TreatmentMetadataItem.setReference("experiment", ExperimentMetadataItem);

                        store(TreatmentMetadataItem);
                    }
                    case "untargeted":
                    {
                        String inhibitorName = treatmentObject.getString("name");
                        String targetGene = treatmentObject.getString("target gene");
                        String concentration = treatmentObject.getString("concentration");     
                        String type = treatmentObject.getString("type");
                        String timePoint = treatmentObject.getString("time point");                            

                        // Save the item
                        Item TreatmentMetadataItem = createItem("NanoporeExperimentTreatment");

                        TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);
                        TreatmentMetadataItem.setAttribute("name", key);

                        if(!StringUtils.isEmpty(concentration) && isDouble(concentration)) {
                            TreatmentMetadataItem.setAttribute("dose_concentration", concentration);
                        }

                        TreatmentMetadataItem.setAttribute("timePoint", timePoint);

                        TreatmentMetadataItem.setReference("experiment", ExperimentMetadataItem);
                        
                        store(TreatmentMetadataItem);
                    }
                }
            } catch (Exception e) {
                LOG.info("Exception in processExperimentTreatments with key: " + key + " - " + e.getMessage());
                continue;
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
