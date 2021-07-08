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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

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
    private Map<String, String> cellLines = new HashMap<String, String>();
    private Map<String, String> RNASeqExperimentMetadatas = new HashMap<String, String>();

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
                    JSONObject experimentJson = (JSONObject)jsonObject.get("experiment");

                    // Get the experiment metadata
                    String experimentName = (String)experimentJson.get("name");
                    String experimentShortName = (String)experimentJson.get("short name");
                    String experimentProject = (String)experimentJson.get("project");
                    String experimentContactPerson = (String)experimentJson.get("contact person");
                    String experimentDate = (String)experimentJson.get("date");
                    String experimentSequencing = (String)experimentJson.get("sequencing");
                    String experimentProvider = (String)experimentJson.get("provider");
                    String experimentDotmaticsReference = (String)experimentJson.get("Dotmatics reference");

                    // Save the item
                    Item ExperimentMetadataItem = createItem("RNASeqExperimentMetadata");                    
                    ExperimentMetadataItem.setAttribute("name", experimentName);
                    ExperimentMetadataItem.setAttribute("shortName", experimentShortName);
                    ExperimentMetadataItem.setAttribute("project", experimentProject);
                    ExperimentMetadataItem.setAttribute("contactPerson", experimentContactPerson);
                    ExperimentMetadataItem.setAttribute("date", experimentDate);
                    ExperimentMetadataItem.setAttribute("sequencing", experimentSequencing);
                    ExperimentMetadataItem.setAttribute("provider", experimentProvider);
                    ExperimentMetadataItem.setAttribute("dotmaticsReference", experimentDotmaticsReference);
                    store(ExperimentMetadataItem);

                    // Process materials
                    JSONObject materialsJson = (JSONObject)jsonObject.get("materials");
                    processRNASeqExperimentMaterials(materialsJson, experimentShortName);

                    // Process treatments
                    JSONObject treatmentsJson = (JSONObject)jsonObject.get("treatments");
                    processRNASeqExperimentTreatments(treatmentsJson, experimentShortName);

                    // Process each comparison individually
                    Map<String, File> filesInDir = readFilesInDir(new File(dataDir.getAbsolutePath() + "/" + experimentShortName));
                    JSONArray experimentComparisons = (JSONArray)jsonObject.get("comparisons");
                    for(int i = 0; i < experimentComparisons.length(); i++) {
                        JSONObject comparison = experimentComparisons.getJSONObject(i);
                        JSONObject treatmentObject = comparison.getJSONObject("treatment");
                        JSONObject controlObject = comparison.getJSONObject("control");
                        String treatmentName = (String)treatmentObject.get("name");
                        String controlName = (String)controlObject.get("name");
                        String fileName = treatmentName + "_vs_" + controlName + "_DESeq2.tsv";

                        if(filesInDir.get(fileName) != null) {
                            File DESeq2File = filesInDir.get(fileName);
                            processRNASeqExperimentDESEQ2(DESeq2File, experimentShortName, treatmentName, controlName);
                        } else {
                            LOG.info("Failed to find DESeq2 file: " + fileName);
                            continue;
                            //throw new RuntimeException("Failed to find DESeq2 file: " + fileName);
                        }
                    }

                    // Process the gene counts
                    String geneCountsFile = "salmon.merged.gene_counts.tsv";
                    if(filesInDir.get(geneCountsFile) != null) {
                        File GeneCountsFile = filesInDir.get(geneCountsFile);
                        processRNASeqExperimentGeneCount(GeneCountsFile, experimentShortName);
                    } else {
                        LOG.info("Failed to find DESeq2 file: " + geneCountsFile);
                        continue;
                        //throw new RuntimeException("Failed to find DESeq2 file: " + geneCountsFile);
                    }

                } catch (JSONException err) {
                    throw new RuntimeException("Failed to read the following JSON: " + line, err);
                }
            }
        }
    }

    private void processRNASeqExperimentMaterials(JSONObject materialsJson, String experimentShortName) {
        // Iterate over each material
        Iterator<String> keys = materialsJson.keys();

        while(keys.hasNext()) {
            String key = keys.next();
            try {
                if (materialsJson.get(key) instanceof JSONObject) {
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
                    JSONObject materialObject = materialsJson.getJSONObject(key).getJSONObject(materialType);

                    switch(materialType) {
                        case "cell line":
                        {
                            String cellLineName = materialObject.getString("name");
                            String cellLineTissue = materialObject.getString("tissue");
                            String cellLineSpecies = materialObject.getString("species");

                            // Save the item
                            Item MaterialMetadataItem = createItem("RNASeqExperimentMaterial");

                            if(!experimentShortName.isEmpty()) {
                                String experimentId = getRNASeqExperimentMetadata(experimentShortName);
                                if (!StringUtils.isEmpty(experimentId)) {   
                                    MaterialMetadataItem.setReference("experiment", experimentId);
                                }
                            }

                            MaterialMetadataItem.setAttribute("materialType", materialType);
                            MaterialMetadataItem.setAttribute("name", cellLineName);
                            MaterialMetadataItem.setAttribute("tissue", cellLineTissue);
                            MaterialMetadataItem.setAttribute("species", cellLineSpecies);
                            store(MaterialMetadataItem);
                        }
                        case "tumour":
                        {
                            String tumourPrimaryDisease = materialObject.getString("primary disease");
                            String tumourDiseaseSubtype = materialObject.getString("disease subtype");
                            String tumourTissue = materialObject.getString("tissue");
                            String tumourSpecies = materialObject.getString("species");

                            // Save the item
                            Item MaterialMetadataItem = createItem("RNASeqExperimentMaterial");

                            if(!experimentShortName.isEmpty()) {
                                String experimentId = getRNASeqExperimentMetadata(experimentShortName);
                                if (!StringUtils.isEmpty(experimentId)) {   
                                    MaterialMetadataItem.setReference("experiment", experimentId);
                                }
                            }

                            MaterialMetadataItem.setAttribute("materialType", materialType);           
                            MaterialMetadataItem.setAttribute("primaryDisease", tumourPrimaryDisease);
                            MaterialMetadataItem.setAttribute("diseaseSubtype", tumourDiseaseSubtype);
                            MaterialMetadataItem.setAttribute("tissue", tumourTissue);
                            MaterialMetadataItem.setAttribute("species", tumourSpecies);
                            store(MaterialMetadataItem);
                        }
                        case "tissue":
                        {
                            String tissueTissue = materialObject.getString("tissue");
                            String tissueSpecies = materialObject.getString("species");

                            // Save the item
                            Item MaterialMetadataItem = createItem("RNASeqExperimentMaterial");

                            if(!experimentShortName.isEmpty()) {
                                String experimentId = getRNASeqExperimentMetadata(experimentShortName);
                                if (!StringUtils.isEmpty(experimentId)) {   
                                    MaterialMetadataItem.setReference("experiment", experimentId);
                                }
                            }

                            MaterialMetadataItem.setAttribute("materialType", materialType);
                            MaterialMetadataItem.setAttribute("tissue", tissueTissue);
                            MaterialMetadataItem.setAttribute("species", tissueSpecies);
                            store(MaterialMetadataItem);
                        }
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
        Iterator<String> keys = treatmentsJson.keys();

        while(keys.hasNext()) {
            String key = keys.next();
            try {
                if (treatmentsJson.get(key) instanceof JSONObject) {
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
                    JSONObject treatmentObject = treatmentsJson.getJSONObject(key).getJSONObject(treatmentType);

                    switch(treatmentType) {
                        case "inhibitor":
                        {
                            String inhibitorName = treatmentObject.getString("name");
                            String targetGene = treatmentObject.getString("target gene");
                            String dotmaticsReference = treatmentObject.getString("Dotmatics reference");
                            String dose = treatmentObject.getString("dose");
                            String timePoint = treatmentObject.getString("time point");                            

                            // Save the item
                            Item TreatmentMetadataItem = createItem("RNASeqExperimentMaterial");

                            if(!experimentShortName.isEmpty()) {
                                String experimentId = getRNASeqExperimentMetadata(experimentShortName);
                                if (!StringUtils.isEmpty(experimentId)) {   
                                    TreatmentMetadataItem.setReference("experiment", experimentId);
                                }
                            }

                            TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);
                            TreatmentMetadataItem.setAttribute("name", key);
                            TreatmentMetadataItem.setAttribute("targetGene", targetGene);
                            TreatmentMetadataItem.setAttribute("dotmaticsReference", dotmaticsReference);

                            if(!StringUtils.isEmpty(dose) && isDouble(dose)) {
                                TreatmentMetadataItem.setAttribute("dose_concentration", dose);
                            }

                            TreatmentMetadataItem.setAttribute("timePoint", timePoint);
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
                            Item TreatmentMetadataItem = createItem("RNASeqExperimentMaterial");

                            if(!experimentShortName.isEmpty()) {
                                String experimentId = getRNASeqExperimentMetadata(experimentShortName);
                                if (!StringUtils.isEmpty(experimentId)) {   
                                    TreatmentMetadataItem.setReference("experiment", experimentId);
                                }
                            }

                            TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);
                            TreatmentMetadataItem.setAttribute("name", key);
                            TreatmentMetadataItem.setAttribute("targetGene", targetGene);

                            if(!StringUtils.isEmpty(concentration) && isDouble(concentration)) {
                                TreatmentMetadataItem.setAttribute("dose_concentration", concentration);
                            }

                            TreatmentMetadataItem.setAttribute("type", type);
                            TreatmentMetadataItem.setAttribute("timePoint", timePoint);
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
                            Item TreatmentMetadataItem = createItem("RNASeqExperimentMaterial");

                            if(!experimentShortName.isEmpty()) {
                                String experimentId = getRNASeqExperimentMetadata(experimentShortName);
                                if (!StringUtils.isEmpty(experimentId)) {   
                                    TreatmentMetadataItem.setReference("experiment", experimentId);
                                }
                            }

                            TreatmentMetadataItem.setAttribute("treatmentType", treatmentType);
                            TreatmentMetadataItem.setAttribute("name", key);

                            if(!StringUtils.isEmpty(concentration) && isDouble(concentration)) {
                                TreatmentMetadataItem.setAttribute("dose_concentration", concentration);
                            }

                            TreatmentMetadataItem.setAttribute("timePoint", timePoint);
                            store(TreatmentMetadataItem);
                        }
                    }
                }
            } catch (Exception e) {
                LOG.info("Exception in processRNASeqExperimentTreatments with key: " + key + " - " + e.getMessage());
                continue;
            }
        }
    }

    private void processRNASeqExperimentDESEQ2(File DESeq2File, String experimentName, String treatmentName, String controlName) throws ObjectStoreException, IOException {
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

                Item IntegratedItem = createItem("RNASeqExperimentDESeq2Result");

                if(!gene.isEmpty()) {
                    String geneId = getGeneId(gene);
                    IntegratedItem.setReference("gene", geneId);
                } else {
                    throw new BuildException("[processRNASeqExperimentDESEQ2] gene was empty: " + experimentName + " - " + fileName);
                }
                               
                IntegratedItem.setAttribute("control", controlName);
                IntegratedItem.setAttribute("treatment", treatmentName);                    

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

                if(!experimentName.isEmpty()) {
                    String experimentId = getRNASeqExperimentMetadata(experimentName);
                    if (!StringUtils.isEmpty(experimentId)) {   
                        IntegratedItem.setReference("experiment", experimentId);
                    }
                }

                store(IntegratedItem);
            }
        }        
    }

    private void processRNASeqExperimentGeneCount(File geneCountsFile, String experimentName) throws ObjectStoreException, IOException {
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

                    if(!StringUtils.isEmpty(count) && isDouble(count)) {
                        IntegratedItem.setAttribute("count", count);
                    }

                    if(!experimentName.isEmpty()) {
                        String experimentId = getRNASeqExperimentMetadata(experimentName);
                        if (!StringUtils.isEmpty(experimentId)) {   
                            IntegratedItem.setReference("experiment", experimentId);
                        }
                    }

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

    public String getRNASeqExperimentMetadata(String identifier) {
        String refId = RNASeqExperimentMetadatas.get(identifier);
        if (refId == null) {
            Item cl = createItem("RNASeqExperimentMetadata");
            cl.setAttribute("shortName", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store RNASeqExperimentMetadata with shortName: " + identifier, e);
            }
            refId = cl.getIdentifier();
            RNASeqExperimentMetadatas.put(identifier, refId);
        }
        return refId;
    }

    public String getCellLine(String identifier) {
        String refId = cellLines.get(identifier);
        if (refId == null) {
            Item cl = createItem("CellLine");
            cl.setAttribute("ShortName", identifier);
            try {
                store(cl);
            } catch (ObjectStoreException e) {
                throw new RuntimeException("failed to store cell line with ShortName: " + identifier, e);
            }
            refId = cl.getIdentifier();
            cellLines.put(identifier, refId);
        }
        return refId;
    }

    private String getGeneId(String primaryIdentifier) throws ObjectStoreException {
        try {
            String resolvedIdentifier = resolveGene(primaryIdentifier);
            String geneId;
            if (StringUtils.isEmpty(resolvedIdentifier)) {
                geneId = genes.get(primaryIdentifier);
                if (geneId == null) {
                    Item gene = createItem("Gene");
                    gene.setAttribute("primaryIdentifier", primaryIdentifier);
                    store(gene);
                    geneId = gene.getIdentifier();
                    genes.put(primaryIdentifier, geneId);
                }
            } else {
                geneId = genes.get(resolvedIdentifier);
                if (geneId == null) {
                    Item gene = createItem("Gene");
                    gene.setAttribute("primaryIdentifier", resolvedIdentifier);
                    store(gene);
                    geneId = gene.getIdentifier();
                    genes.put(resolvedIdentifier, geneId);
                }
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
