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

import org.intermine.dataconversion.MockItemWriter;
import org.intermine.model.fulldata.Item;
import org.intermine.dataconversion.ItemsTestCase;
import org.intermine.metadata.Model;

import java.io.File;
import java.util.HashMap;

public class DepmapSampleInfoConverterTest extends ItemsTestCase
{
    Model model = Model.getInstanceByName("genomic");
    DepmapSampleInfoConverter converter;
    MockItemWriter itemWriter;

    public DepmapSampleInfoConverterTest(String arg) {
        super(arg);
    }

    public void setUp() throws Exception {
        super.setUp();
        itemWriter = new MockItemWriter(new HashMap<String, Item>());
        converter = new DepmapSampleInfoConverter(itemWriter, model);
    }

    public void testProcess() throws Exception {
        File tmp = new File(getClass().getClassLoader().getResource("sample_info.csv").toURI());
        File dataDirectory = tmp.getParentFile();

        System.out.println(dataDirectory.getAbsolutePath());

        converter.process(dataDirectory);
        converter.close();

        writeItemsFile(itemWriter.getItems(), "depmap-sample-info-items.xml");
    }
}
