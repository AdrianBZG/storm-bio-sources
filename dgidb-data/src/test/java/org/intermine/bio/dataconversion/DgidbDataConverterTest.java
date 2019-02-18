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
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.HashMap;
import java.util.Set;

public class DgidbDataConverterTest extends ItemsTestCase
{
    Model model = Model.getInstanceByName("genomic");

    DgidbDataConverter converter;
    MockItemWriter itemWriter;

    public DgidbDataConverterTest(String arg) {
        super(arg);
    }

    public void setUp() throws Exception {
        super.setUp();
        itemWriter = new MockItemWriter(new HashMap<String, Item>());
        converter = new DgidbDataConverter(itemWriter, model);
    }


    public void testProcess() throws Exception {
        File tmp = new File(getClass().getClassLoader().getResource("drugs.tsv").toURI());
        File dataDirectory = tmp.getParentFile();

        System.out.println(dataDirectory.getAbsolutePath());

        converter.process(dataDirectory);
        converter.close();
        // uncomment to write out a new target items file
        writeItemsFile(itemWriter.getItems(), "dgidb-items.xml");

        Set<org.intermine.xml.full.Item> expected = readItemSet("DgidbDataConverterTest-tgt.xml");
        assertEquals(expected, itemWriter.getItems());
    }
}
