package infodynamics.measures.discrete;

import static org.junit.Assert.assertThrows;

import org.junit.Test;

import infodynamics.measures.discrete.InfoMeasureCalculatorDiscrete.State;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.List;

public class EntropyTester extends TestCase {
    
    private EntropyCalculatorDiscrete calc;
    private static double ERR_ALLOWANCE = 0.02;

    public void setup() {
        calc = new EntropyCalculatorDiscrete();
    }

    @Test
    public void testHashTableUseBasic() throws Exception {
        setup();

        // 1. check the state machine is in the correct state
        assertEquals(State.SETTING_PROPERTIES, calc.currentState);

        calc.initialise();
        String[] obs = new String[] {"a", "b", "b"};
        calc.addObservations(obs);

        // 2. check that the integer range is infact unknown
        // 3. check the state machine
        assertEquals(false, calc.knownIntegerRange);
        assertEquals(State.ADDING_OBSERVATIONS, calc.currentState);
        
        // 4. check the hash table and it's values
        assertEquals((Integer) 1, calc.hashedStateCount.get("a"));
        assertEquals((Integer) 2, calc.hashedStateCount.get("b"));

        double ent = calc.computeAverageLocalOfObservations();

        // 5. check result of entropy calculations - allowing small error margin
        assertTrue(Math.abs(ent-0.918)/0.918 <= ERR_ALLOWANCE);
        assertEquals(State.COMPUTING, calc.currentState);
    }

    @Test
    public void testArrayUseBasic() throws Exception {
        setup();
        calc.setProperty("ALPHABET_SIZE", "2");
        
        assertEquals(State.SETTING_PROPERTIES, calc.currentState);
        assertEquals(true, calc.knownIntegerRange);

        int[] obs = new int[] {0,0,1,1};
        calc.addObservations(obs);

        assertEquals(State.ADDING_OBSERVATIONS, calc.currentState);

        double ent = calc.computeAverageLocalOfObservations();
        assertEquals(1.0, ent);
        assertEquals(State.COMPUTING, calc.currentState);
    }

    @Test
    public void testArrayUseOverflowAlphaSize() throws Exception {
        setup();
        calc.setProperty("ALPHABET_SIZE", "2");

        int[] obs = new int[] { 0, 1, 2 };
        assertThrows(RuntimeException.class, () -> {
            calc.addObservations(obs);
        });
    }

    @Test
    public void testArrayUseOverflowAlphaSize2() throws Exception {
        setup();
        calc.setProperty("ALPHABET_SIZE", "2");

        Integer[] obs = new Integer[] { 0, 1, 2 };
        assertThrows(RuntimeException.class, () -> {
            calc.addObservations(obs);
        });
    }

    @Test
    public void testStringIntegerInterpretation() throws Exception {
        setup();
        calc.setProperty("ALPHABET_SIZE", "2");
        
        String[] obs = new String[] { "0", "1", "1" };
        calc.addObservations(obs);

        double ent = calc.computeAverageLocalOfObservations();

        // ensure we did use the array implementation
        assertEquals(calc.knownIntegerRange, true);
        assertNotNull(calc.stateCount);
        assertTrue(calc.hashedStateCount.isEmpty());

        // ensure it gets correct output - allowing small error margin
        assertTrue(Math.abs(ent-0.918)/0.918 <= ERR_ALLOWANCE);
    }

    @Test
    public void testUninterpretableInput() throws Exception {
        setup();
        calc.setProperty("ALPHABET_SIZE", "2");
        
        String[] obs = new String[] { "0", "1", "a" };

        // we specified an alphabet size, but cannot interpret "a" as an integer
        assertThrows(RuntimeException.class, () ->{
            calc.addObservations(obs);
        });
    }

    @Test
    public void testNegativeNums() throws Exception {
        setup();

        calc.setProperty("ALPHABET_SIZE", "4");
        int[] obs = new int[] { 0, 1, -1, 2 };

        // Negative numbers are invalid observations.
        assertThrows(RuntimeException.class, () -> {
            calc.addObservations(obs);
        });
    }

    @Test
    public void testExceedMaxAlphaSize() throws Exception {
        setup();

        // TODO -- This test is under the assumption of the temporary max alphabet size of 100
        // when this number is solidified, this test needs to change (it will fail when it does change.)
        int[] obs = new int[] { 99, 100, 101 };
        calc.addObservations(obs);

        assertFalse(calc.hashedStateCount.isEmpty());
        assertNull(calc.stateCount);
    }

    @Test
    public void testMultiDimensionalIntObservations() throws Exception {
        setup();

        calc.setProperty("NUM_DIMENSIONS", "3");

        int[][] obs = new int[][] {
            {1, 2, 3},
            {4, 5, 6},
            {1, 2, 3}
        };

        calc.addObservations(obs);

        // TODO -- compute combined values


        // assertEquals((Integer) 2, calc.hashedStateCount.get(test1));
        // assertEquals((Integer) 1, calc.hashedStateCount.get(test2));
    }

    @Test
    public void testMultiDimensionalObjObservations() throws Exception {
        setup();

        calc.setProperty("NUM_DIMENSIONS", "3");

        Object[][] obs = new Object[][] {
            {"A", "B", "C"},
            {"D", "E", "F"},
            {"A", "B", "C"}
        };

        calc.addObservations(obs);

        List<Object> test1 = Arrays.asList(obs[0]);
        List<Object> test2 = Arrays.asList(obs[1]);

        assertEquals((Integer) 2, calc.hashedStateCount.get(test1));
        assertEquals((Integer) 1, calc.hashedStateCount.get(test2));
    }

    // TODO -- implement this test.
    @Test
    public void testZeroDimensions() throws Exception {
        setup();

        calc.setProperty("NUM_DIMENSIONS", "0");
    }
}
