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
public class StormRnaseqDataConverter extends BioDirectoryConverter
{
    private static final String DATASET_TITLE = "STORM RNA-Seq Data";
    private static final String DATA_SOURCE_NAME = "STORM RNA-Seq Data";

    private static final String TAXON_ID = "9606"; // Human Taxon ID

    private Map<String, String> genes = new HashMap<String, String>();
    private Map<String, String> resolvedGenes = new HashMap<String, String>();
    private Map<String, String> unresolvableGenes = new HashMap<String, String>();
    private Map<String, Item> experiments = new HashMap<>();
    private Map<String, Item> materials = new HashMap<>();
    private Map<String, Item> treatments = new HashMap<>();
    private Map<String, Item> conditions = new HashMap<>();

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
        //    
    }

    private Map<String, Integer> getColumnIndexes(String[] header, String fileType) {
        Map<String, Integer> indexes = new HashMap<String, Integer>();

        if(fileType == "DESEQ2") {
            //20: Ensembl	Entrez	Gene	Description	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
            //19: ensembl	entrez	symbol	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
            for(int i = 0; i < header.length; i++) {
                String column = header[i];
                if(column.toLowerCase().contains("ensembl")) {
                    indexes.put("gene", i);
                    continue;
                } else if(column.toLowerCase().contains("basemean")) {
                    indexes.put("baseMean", i);
                    continue;
                } else if(column.toLowerCase().contains("log2foldchange")) {
                    indexes.put("log2FoldChange", i);
                    continue;
                } else if(column.toLowerCase().contains("lfcse")) {
                    indexes.put("lfcSE", i);
                    continue;
                } else if(column.toLowerCase().contains("stat")) {
                    indexes.put("stat", i);
                    continue;
                } else if(column.toLowerCase().contains("pvalue")) {
                    indexes.put("pvalue", i);
                    continue;
                } else if(column.toLowerCase().contains("padj")) {
                    indexes.put("padj", i);
                    continue;
                }
            }
        } else if(fileType == "SampleInfo") {
            //20: Run,Sample,Condition,Celline,IFN_gamma,Compound,Concentration µM,Replicate,Timepoint,Lexogen_ID,Filename
            //19: Run,Condition,Celline,IFN gamma,Compound,Concentration,Replicate,Timepoint,Sample number,Sample,File name
            for(int i = 0; i < header.length; i++) {
                String column = header[i];
                if(column.toLowerCase().contains("run")) {
                    indexes.put("run", i);
                    continue;
                } else if(column.toLowerCase().equals("sample")) {
                    indexes.put("sample", i);
                    continue;
                } else if(column.toLowerCase().contains("condition") || column.toLowerCase().contains("treatment")) {
                    indexes.put("treatment", i);
                    continue;
                } else if(column.toLowerCase().contains("cellline") || column.toLowerCase().contains("celline")) {
                    indexes.put("cellLine", i);
                    continue;
                } else if(column.toLowerCase().contains("ifn_gamma") || column.toLowerCase().contains("ifn gamma")) {
                    indexes.put("IFN_gamma", i);
                    continue;
                } else if(column.toLowerCase().contains("compound")) {
                    indexes.put("compound", i);
                    continue;
                } else if(column.toLowerCase().contains("concentration")) {
                    indexes.put("concentration", i);
                    continue;
                } else if(column.toLowerCase().contains("replicate")) {
                    indexes.put("replicate", i);
                    continue;
                } else if(column.toLowerCase().contains("timepoint")) {
                    indexes.put("timepoint", i);
                    continue;
                }
            }
        } else {
            throw new BuildException("Unknown fileType for getColumnIndexes: " + fileType);
        }

        return indexes;
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

                    LOG.info("StormRnaseqDataConverter [processConfigFile] - Processing 1: " + experimentShortName);

                    // Save the item
                    Item ExperimentMetadataItem = createItem("RNASeqExperimentMetadata");

                    if(!experimentName.isEmpty()) {                 
                        ExperimentMetadataItem.setAttribute("name", experimentName);
                    }

                    if(!experimentShortName.isEmpty()) {
                        ExperimentMetadataItem.setAttribute("shortName", experimentShortName);
                    }

                    if(!experimentProject.isEmpty()) {
                        ExperimentMetadataItem.setAttribute("project", experimentProject);
                    }

                    if(!experimentContactPerson.isEmpty()) {
                        ExperimentMetadataItem.setAttribute("contactPerson", experimentContactPerson);
                    }

                    if(!experimentDate.isEmpty()) {
                        ExperimentMetadataItem.setAttribute("date", experimentDate);
                    }

                    if(!experimentSequencing.isEmpty()) {
                        ExperimentMetadataItem.setAttribute("sequencing", experimentSequencing);
                    }

                    if(!experimentProvider.isEmpty()) {
                        ExperimentMetadataItem.setAttribute("provider", experimentProvider);
                    }

                    if(!experimentDotmaticsReference.isEmpty()) {
                        ExperimentMetadataItem.setAttribute("dotmaticsReference", experimentDotmaticsReference);
                    }
                    //store(ExperimentMetadataItem);

                    String experimentKey = experimentShortName;
                    if (!experiments.containsKey(experimentKey)) {
                        experiments.put(experimentKey, ExperimentMetadataItem);
                    }

                    // Process materials
                    LOG.info("StormRnaseqDataConverter [processConfigFile] - Processing 2: " + experimentShortName);
                    JSONObject materialsJson = jsonObject.getJSONObject("materials");
                    processRNASeqExperimentMaterials(materialsJson, experimentShortName);

                    // Process treatments
                    LOG.info("StormRnaseqDataConverter [processConfigFile] - Processing 3: " + experimentShortName);
                    JSONObject treatmentsJson = jsonObject.getJSONObject("treatments");
                    processRNASeqExperimentTreatments(treatmentsJson, experimentShortName);

                    // Process conditions
                    JSONObject conditionsJson = jsonObject.getJSONObject("conditions");
                    processRNASeqExperimentConditions(conditionsJson, experimentShortName);

                    // Process each comparison individually
                    LOG.info("StormRnaseqDataConverter [processConfigFile] - Processing 4: " + experimentShortName);
                    Map<String, File> filesInDir = readFilesInDir(new File(dataDir.getAbsolutePath() + "/" + experimentShortName));
                    JSONArray experimentComparisons = (JSONArray)jsonObject.get("comparisons");
                    for(int i = 0; i < experimentComparisons.length(); i++) {
                        JSONObject comparison = experimentComparisons.getJSONObject(i);
                        JSONObject treatmentObject = comparison.getJSONObject("treatment");
                        JSONObject controlObject = comparison.getJSONObject("control");
                        String treatmentName = treatmentObject.getString("name");
                        String controlName = controlObject.getString("name");
                        String fileName = treatmentName + "_vs_" + controlName + "_DESeq2.tsv";

                        if(filesInDir.get(fileName) != null) {
                            File DESeq2File = filesInDir.get(fileName);
                            processRNASeqExperimentComparison(DESeq2File, experimentShortName, treatmentName, controlName);
                        } else {
                            LOG.info("Failed to find DESeq2 file: " + fileName);
                            continue;
                            //throw new RuntimeException("Failed to find DESeq2 file: " + fileName);
                        }
                    }

                    // Process the gene counts
                    LOG.info("StormRnaseqDataConverter [processConfigFile] - Processing 5: " + experimentShortName);
                    String geneCountsFile = "salmon.merged.gene_counts.tsv";
                    if(filesInDir.get(geneCountsFile) != null) {
                        File GeneCountsFile = filesInDir.get(geneCountsFile);
                        processRNASeqExperimentGeneCount(GeneCountsFile, experimentShortName);
                    } else {
                        LOG.info("Failed to find DESeq2 file: " + geneCountsFile);
                        continue;
                        //throw new RuntimeException("Failed to find DESeq2 file: " + geneCountsFile);
                    }

                    // Store the experiment
                    LOG.info("StormRnaseqDataConverter [processConfigFile] - Processing 6: " + experimentShortName);
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
            throw new RuntimeException("Error storing StormRNASeqExperiment ", e);
        }
    }

    private void processRNASeqExperimentMaterials(JSONObject materialsJson, String experimentShortName) {
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
                
                if(materialTypeKeysJSON.has("cell line")) {
                    String materialType = "cell line";
                
                    // Get the material object
                    JSONObject materialObject = materialTypeKeysJSON.getJSONObject(materialType);

                    String cellLineName = "";
                    if(materialObject.has("name")) {
                        cellLineName = materialObject.getString("name");
                    }

                    String cellLineTissue = "";
                    if(materialObject.has("tissue")) {
                        cellLineTissue = materialObject.getString("tissue");
                    }

                    String cellLineSpecies = "";
                    if(materialObject.has("species")) {
                        cellLineSpecies = materialObject.getString("species");
                    }

                    // Save the item
                    Item MaterialMetadataItem = createItem("RNASeqExperimentMaterial");                       

                    if(!materialType.isEmpty()) {
                        MaterialMetadataItem.setAttribute("materialType", materialType);
                    }

                    if(!cellLineName.isEmpty()) {
                        MaterialMetadataItem.setAttribute("name", cellLineName);
                    }

                    if(!cellLineTissue.isEmpty()) {
                        MaterialMetadataItem.setAttribute("tissue", cellLineTissue);
                    }

                    if(!cellLineSpecies.isEmpty()) {
                        MaterialMetadataItem.setAttribute("species", cellLineSpecies);
                    }

                    MaterialMetadataItem.setReference("experiment", ExperimentMetadataItem);

                    store(MaterialMetadataItem);

                    if (!materials.containsKey(materialName)) {
                        materials.put(materialName, MaterialMetadataItem);
                    }
                } else if(materialTypeKeysJSON.has("tumour")) {
                    String materialType = "tumour";

                    // Get the material object
                    JSONObject materialObject = materialTypeKeysJSON.getJSONObject(materialType);

                    String tumourPrimaryDisease = "";
                    if(materialObject.has("primary disease")) {
                        tumourPrimaryDisease = materialObject.getString("primary disease");
                    }

                    String tumourDiseaseSubtype = "";
                    if(materialObject.has("disease subtype")) {
                        tumourDiseaseSubtype = materialObject.getString("disease subtype");
                    }

                    String tumourTissue = "";
                    if(materialObject.has("tissue")) {
                        tumourTissue = materialObject.getString("tissue");
                    }

                    String tumourSpecies = "";
                    if(materialObject.has("species")) {
                        tumourSpecies = materialObject.getString("species");
                    }

                    // Save the item
                    Item MaterialMetadataItem = createItem("RNASeqExperimentMaterial");

                    if(!materialType.isEmpty()) {
                        MaterialMetadataItem.setAttribute("materialType", materialType);
                    }
                    
                    if(!tumourPrimaryDisease.isEmpty()) {
                        MaterialMetadataItem.setAttribute("primaryDisease", tumourPrimaryDisease);
                    }

                    if(!tumourDiseaseSubtype.isEmpty()) {
                        MaterialMetadataItem.setAttribute("diseaseSubtype", tumourDiseaseSubtype);
                    }

                    if(!tumourTissue.isEmpty()) {
                        MaterialMetadataItem.setAttribute("tissue", tumourTissue);
                    }

                    if(!tumourSpecies.isEmpty()) {
                        MaterialMetadataItem.setAttribute("species", tumourSpecies);
                    }

                    MaterialMetadataItem.setReference("experiment", ExperimentMetadataItem);

                    store(MaterialMetadataItem);

                    if (!materials.containsKey(materialName)) {
                        materials.put(materialName, MaterialMetadataItem);
                    }
                } else if(materialTypeKeysJSON.has("tissue")) {
                    String materialType = "tissue";

                    // Get the material object
                    JSONObject materialObject = materialTypeKeysJSON.getJSONObject(materialType);

                    String tissueTissue = "";
                    if(materialObject.has("tissue")) {
                        tissueTissue = materialObject.getString("tissue");
                    }

                    String tissueSpecies = "";
                    if(materialObject.has("species")) {
                        tissueSpecies = materialObject.getString("species");
                    }

                    // Save the item
                    Item MaterialMetadataItem = createItem("RNASeqExperimentMaterial");

                    if(!materialType.isEmpty()) {
                        MaterialMetadataItem.setAttribute("materialType", materialType);
                    }

                    if(!tissueTissue.isEmpty()) {
                        MaterialMetadataItem.setAttribute("tissue", tissueTissue);
                    }

                    if(!tissueSpecies.isEmpty()) {
                        MaterialMetadataItem.setAttribute("species", tissueSpecies);
                    }

                    MaterialMetadataItem.setReference("experiment", ExperimentMetadataItem);   

                    store(MaterialMetadataItem);

                    if (!materials.containsKey(materialName)) {
                        materials.put(materialName, MaterialMetadataItem);
                    }
                }
                              
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentMaterials with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processRNASeqExperimentTreatments(JSONObject treatmentsJson, String experimentShortName) {
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

                if(treatmentTypeKeysJSON.has("inhibitor")) {
                    String treatmentType = "inhibitor";
                
                    // Get the material object
                    JSONObject treatmentObject = treatmentTypeKeysJSON.getJSONObject(treatmentType);

                    String inhibitorName = "";
                    if(treatmentObject.has("name")) {
                        inhibitorName = treatmentObject.getString("name");
                    }

                    String targetGene = "";
                    if(treatmentObject.has("target gene")) {
                        targetGene = treatmentObject.getString("target gene");
                    }

                    String dotmaticsReference = "";
                    if(treatmentObject.has("Dotmatics reference")) {
                        dotmaticsReference = treatmentObject.getString("Dotmatics reference");
                    }

                    String dose = "";
                    if(treatmentObject.has("dose")) {
                        dose = treatmentObject.getString("dose");
                    }

                    String timePoint = "";
                    if(treatmentObject.has("time point")) {
                        timePoint = treatmentObject.getString("time point");
                    }                     

                    // Save the item
                    Item TreatmentMetadataItem = createItem("RNASeqExperimentTreatment");

                    TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);                    

                    if(!key.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("name", key);
                    }

                    if(!targetGene.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("targetGene", targetGene);
                    }

                    if(!dotmaticsReference.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("dotmaticsReference", dotmaticsReference);
                    }

                    if(!StringUtils.isEmpty(dose) && isDouble(dose)) {
                        TreatmentMetadataItem.setAttribute("dose_concentration", dose);
                    }

                    if(!timePoint.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("timePoint", timePoint);
                    }

                    TreatmentMetadataItem.setReference("experiment", ExperimentMetadataItem);

                    store(TreatmentMetadataItem);

                    if (!treatments.containsKey(treatmentName)) {
                        treatments.put(treatmentName, TreatmentMetadataItem);
                    }
                    
                } else if(treatmentTypeKeysJSON.has("knock-down")) {
                    String treatmentType = "knock-down";                
                    // Get the material object
                    JSONObject treatmentObject = treatmentTypeKeysJSON.getJSONObject(treatmentType);

                    String inhibitorName = "";
                    if(treatmentObject.has("name")) {
                        inhibitorName = treatmentObject.getString("name");
                    }

                    String targetGene = "";
                    if(treatmentObject.has("target gene")) {
                        targetGene = treatmentObject.getString("target gene");
                    }

                    String concentration = "";
                    if(treatmentObject.has("concentration")) {
                        concentration = treatmentObject.getString("concentration");
                    }

                    String type = "";
                    if(treatmentObject.has("type")) {
                        type = treatmentObject.getString("type");
                    }

                    String timePoint = "";
                    if(treatmentObject.has("time point")) {
                        timePoint = treatmentObject.getString("time point");
                    }                     

                    // Save the item
                    Item TreatmentMetadataItem = createItem("RNASeqExperimentTreatment");
            
                    TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);                    

                    if(!key.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("name", key);
                    }

                    if(!targetGene.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("targetGene", targetGene);
                    }

                    if(!StringUtils.isEmpty(concentration) && isDouble(concentration)) {
                        TreatmentMetadataItem.setAttribute("dose_concentration", concentration);
                    }

                    if(!type.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("type", type);
                    }

                    if(!timePoint.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("timePoint", timePoint);
                    }

                    TreatmentMetadataItem.setReference("experiment", ExperimentMetadataItem);

                    store(TreatmentMetadataItem);

                    if (!treatments.containsKey(treatmentName)) {
                        treatments.put(treatmentName, TreatmentMetadataItem);
                    }
                }
                else if(treatmentTypeKeysJSON.has("untargeted")) {
                    String treatmentType = "untargeted";                
                    // Get the material object
                    JSONObject treatmentObject = treatmentTypeKeysJSON.getJSONObject(treatmentType);

                    String inhibitorName = "";
                    if(treatmentObject.has("name")) {
                        inhibitorName = treatmentObject.getString("name");
                    }

                    String targetGene = "";
                    if(treatmentObject.has("target gene")) {
                        targetGene = treatmentObject.getString("target gene");
                    }

                    String concentration = "";
                    if(treatmentObject.has("concentration")) {
                        concentration = treatmentObject.getString("concentration");
                    }
                    
                    String type = "";
                    if(treatmentObject.has("type")) {
                        type = treatmentObject.getString("type");
                    }

                    String timePoint = "";
                    if(treatmentObject.has("time point")) {
                        timePoint = treatmentObject.getString("time point");
                    }                     

                    // Save the item
                    Item TreatmentMetadataItem = createItem("RNASeqExperimentTreatment");
      
                    TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);

                    if(!key.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("name", key);
                    }

                    if(!StringUtils.isEmpty(concentration) && isDouble(concentration)) {
                        TreatmentMetadataItem.setAttribute("dose_concentration", concentration);
                    }

                    if(!timePoint.isEmpty()) {
                        TreatmentMetadataItem.setAttribute("timePoint", timePoint);
                    }

                    TreatmentMetadataItem.setReference("experiment", ExperimentMetadataItem);

                    store(TreatmentMetadataItem);

                    if (!treatments.containsKey(treatmentName)) {
                        treatments.put(treatmentName, TreatmentMetadataItem);
                    }
                }                
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentTreatments with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processRNASeqExperimentConditions(JSONObject conditionsJson, String experimentShortName) {
        // Iterate over each condition
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        Iterator<String> keys = conditionsJson.keys();

        while(keys.hasNext()) {
            String key = keys.next();
            try {
                // Now we have the condition
                String conditionName = key;

                // There should only be one key under this
                JSONObject conditionsKeysJSON = conditionsJson.getJSONObject(conditionName);

                String materialName = conditionsKeysJSON.getString("material");

                // Samples
                ArrayList<String> samplesArray = new ArrayList<String>();
                JSONObject samplesJson = conditionsKeysJSON.getJSONObject("samples");
                Iterator<String> samplesKeys = samplesJson.keys();
                while(samplesKeys.hasNext()) {
                    String sampleKey = samplesKeys.next();
                    samplesArray.add(sampleKey);
                }

                String samples = String.join(", ", samplesArray);

                // Treatments
                ArrayList<String> treatmentsArray = new ArrayList<String>();
                JSONArray treatmentsJson = conditionsKeysJSON.getJSONArray("treatments");
                for(int i = 0; i < treatmentsJson.length(); i++) {
                    treatmentsArray.add(treatmentsJson.getString(i));
                }

                String treatments = String.join(", ", treatmentsArray);

                // Save the item
                Item ConditionMetadataItem = createItem("RNASeqExperimentCondition");

                if(!conditionName.isEmpty()) {
                    ConditionMetadataItem.setAttribute("name", conditionName);
                } else {
                    continue;
                }

                if(!treatments.isEmpty()) {
                    ConditionMetadataItem.setAttribute("treatments", treatments);
                }

                if(!samples.isEmpty()) {
                    ConditionMetadataItem.setAttribute("samples", samples);
                }

                if(!materialName.isEmpty()) {
                    if (!materials.containsKey(materialName)) {
                        ConditionMetadataItem.setReference("material", materials.get(materialName));
                    }
                }

                ConditionMetadataItem.setReference("experiment", ExperimentMetadataItem);

                store(ConditionMetadataItem);

                if (!conditions.containsKey(conditionName)) {
                    conditions.put(conditionName, ConditionMetadataItem);
                }

    
                
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentConditions with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processRNASeqExperimentComparison(File DESeq2File, String experimentShortName, String treatmentName, String controlName) throws ObjectStoreException, IOException {
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        String fileName = DESeq2File.getName();

        if(fileName.endsWith("_DESeq2.tsv")) {
            String fileAbsPath = DESeq2File.getAbsolutePath();

            Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
            String[] firstLine = (String[]) lineIter.next();
            Map<String, Integer> columnIndexes = getColumnIndexes(firstLine, "DESEQ2");

            while (lineIter.hasNext()) {
                String[] line = (String[]) lineIter.next();

                String gene = line[columnIndexes.get("gene").intValue()].split("\\.")[0];
                String baseMean = line[columnIndexes.get("baseMean").intValue()];
                String log2FoldChange = line[columnIndexes.get("log2FoldChange").intValue()];
                String lfcSE = line[columnIndexes.get("lfcSE").intValue()];
                String stat = line[columnIndexes.get("stat").intValue()];
                String pvalue = line[columnIndexes.get("pvalue").intValue()];
                String padj = line[columnIndexes.get("padj").intValue()];

                Item IntegratedItem = createItem("RNASeqExperimentComparison");

                if(!gene.isEmpty()) {
                    if(unresolvableGenes.get(gene) != null) {
                        continue;
                    }
                    String geneId = getGeneId(gene);
                    if(geneId == null) {
                        continue;
                    }
                    IntegratedItem.setReference("gene", geneId);
                } else {
                    throw new BuildException("[processRNASeqExperimentComparison] gene was empty: " + fileName);
                }
                               
                if(conditions.containsKey(controlName)) {
                    IntegratedItem.setReference("control", conditions.get(controlName));
                } else {
                    continue;
                }

                if(conditions.containsKey(treatmentName)) {
                    IntegratedItem.setReference("treatment", conditions.get(treatmentName));
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

                IntegratedItem.setReference("experiment", ExperimentMetadataItem);         

                store(IntegratedItem);
            }
        }
    }

    private void processRNASeqExperimentGeneCount(File geneCountsFile, String experimentShortName) throws ObjectStoreException, IOException {
        Item ExperimentMetadataItem = experiments.get(experimentShortName);
        String fileAbsPath = geneCountsFile.getAbsolutePath();
        Iterator<?> lineIter = FormattedTextParser.parseTabDelimitedReader(new FileReader(fileAbsPath));
        String[] firstLine = (String[]) lineIter.next();

        ArrayList<String> runs = new ArrayList<String>();
        for(int i = 2; i < firstLine.length; i++) {
            String run = firstLine[i];
            runs.add(run);
        }

        while (lineIter.hasNext()) {

            String[] line = (String[]) lineIter.next();
            String gene = line[1];
            try {
                for(int i = 2; i < line.length; i++) {
                    String count = line[i];
                    String runForThisItem = runs.get(i-2);
                    Item IntegratedItem = createItem("RNASeqExperimentGeneCount");
                    if(!gene.isEmpty()) {
                        if(unresolvableGenes.get(gene) != null) {
                            continue;
                        }
                        String geneId = getGeneId(gene);
                        if (geneId == null) {
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

                    if(!StringUtils.isEmpty(count) && isDouble(count)) {
                        IntegratedItem.setAttribute("count", count);
                    }

                    IntegratedItem.setReference("experiment", ExperimentMetadataItem);

                    store(IntegratedItem);
                }
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentGeneCount with gene: " + gene + " - " + e.getMessage());
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

    private boolean isDouble(String str) {
        try {
            double x = Double.parseDouble(str);
            //if (x == (int) x)
            //    return false;
            return true;
        }
        catch(NumberFormatException e) {
            return false;
        }
    }
}
