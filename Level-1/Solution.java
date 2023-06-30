import java.util.HashMap;
import java.util.HashSet;

public class Solution {
    public static int solution(int[] x, int[] y) {
        /* Explanation:
            We know that x and y are similar, but not the same. 
            One list will have ONE extra value in it.
            We can loop through both lists and add the values to a hashmap
            We can set the key to the ID and the value to the occurances
            We can then loop through the hashmap and find the key with the value of 1
            We can then return that key.
            But, since that requires two loops, 
            What if we just removed all elements when they are visited twice?
        */
        // Check if x is larger than y
        // This is used later to fix an edge case in the for loop
        boolean xIsLarger = x.length > y.length;
        // Create a map for the IDs, with the key as the ID 
        // and the value as the occurance count
        HashMap<Integer, Integer> IDs = new HashMap<Integer, Integer>();
        // Create a set for the IDs that have been visited twice,
        // which we KNOW cannot be the extra
        HashSet<Integer> blacklistedIDs = new HashSet<Integer>();
        
        // Loop through the lists and add the values to the hashmap
        int maxIdx = xIsLarger ? x.length : y.length;
        for (int i = 0; i < (maxIdx); i++) 
        {
            // Check if we have reached the end of the smaller list
            // only add the last value to the hashmap
            if (i == (maxIdx)-1)
            {
                // set "last value" to the value of the larger list
                // to remove the requirement for redundant code
                int lastValue = xIsLarger ? x[i] : y[i];
                
                safeAdd(lastValue, IDs, blacklistedIDs);
                safeRemove(lastValue, IDs, blacklistedIDs);
                
                break;
            }
            // i is in range of both lists,
            // add them both to the hashmap
            // if they are not already whitelisted
            else
            {
                safeAdd(x[i], IDs, blacklistedIDs);
                safeAdd(y[i], IDs, blacklistedIDs);
            }
            // Check if the value has been visited twice
            // This means that it is not the extra value, 
            // and we can blacklist it.

            // It is important to note that the extra value will 
            // ALWAYS be unique, and thus will never be visited twice
            // i.e. if it is in the same list twice.
            safeRemove(x[i], IDs, blacklistedIDs);
            safeRemove(y[i], IDs, blacklistedIDs);
        }
        return IDs.keySet().iterator().next();
    }

    /**
     * Removes the value from the IDs hashmap if it has been visited twice
     * @param value          the value to check against the blacklist
     * @param IDs            the hashmap of IDs
     * @param blacklistedIDs the set of blacklisted IDs
     */
    private static void safeRemove(int value, HashMap<Integer, Integer> IDs, HashSet<Integer> blacklistedIDs) {
        if (IDs.getOrDefault(value, 0) == 2)
        {
            blacklistedIDs.add(value);
            IDs.remove(value);
        }
    }

    /**
     * Adds the value to the IDs hashmap if it is not blacklisted
     * @param value          the value to add
     * @param IDs            the hashmap of IDs
     * @param blacklistedIDs the set of blacklisted IDs
     */
    private static void safeAdd(int value, HashMap<Integer, Integer> IDs, HashSet<Integer> blacklistedIDs) {
        if (!blacklistedIDs.contains(value))
        {
            IDs.put(value, IDs.getOrDefault(value, 0) + 1);
        }
    }
}