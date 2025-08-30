import org.json.simple.*;
import org.json.simple.parser.*;
import java.io.*;
import java.math.BigInteger;
import java.util.*;

public class secrete {

    // Exact rational number using BigInteger
    static class Rational {
        BigInteger n, d; // numerator, denominator (>0)
        Rational(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("zero denom");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            if (!g.equals(BigInteger.ONE)) { n = n.divide(g); d = d.divide(g); }
            this.n = n; this.d = d;
        }
        Rational(BigInteger v) { this(v, BigInteger.ONE); }
        Rational(long v) { this(BigInteger.valueOf(v), BigInteger.ONE); }
        static Rational zero() { return new Rational(BigInteger.ZERO, BigInteger.ONE); }
        static Rational one() { return new Rational(BigInteger.ONE, BigInteger.ONE); }
        boolean isZero() { return n.equals(BigInteger.ZERO); }

        Rational add(Rational o) {
            BigInteger nn = n.multiply(o.d).add(o.n.multiply(d));
            BigInteger dd = d.multiply(o.d);
            return new Rational(nn, dd);
        }
        Rational sub(Rational o) {
            BigInteger nn = n.multiply(o.d).subtract(o.n.multiply(d));
            BigInteger dd = d.multiply(o.d);
            return new Rational(nn, dd);
        }
        Rational mul(Rational o) {
            return new Rational(n.multiply(o.n), d.multiply(o.d));
        }
        Rational div(Rational o) {
            if (o.n.equals(BigInteger.ZERO)) throw new ArithmeticException("div by zero");
            return new Rational(n.multiply(o.d), d.multiply(o.n));
        }
        Rational neg() { return new Rational(n.negate(), d); }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Rational)) return false;
            Rational r = (Rational) o;
            return n.equals(r.n) && d.equals(r.d);
        }
        @Override
        public int hashCode(){ return Objects.hash(n,d); }

        @Override
        public String toString() {
            if (d.equals(BigInteger.ONE)) return n.toString();
            return n.toString() + "/" + d.toString();
        }
    }

    // Gauss-Jordan solve augmented matrix (k x (k+1)). Returns solution array or null if singular.
    static Rational[] solveAugmented(Rational[][] aug) {
        int n = aug.length;
        int m = n + 1;
        for (int col = 0; col < n; col++) {
            // find pivot
            int piv = -1;
            for (int r = col; r < n; r++) {
                if (!aug[r][col].isZero()) { piv = r; break; }
            }
            if (piv == -1) return null; // singular
            if (piv != col) {
                Rational[] tmp = aug[col]; aug[col] = aug[piv]; aug[piv] = tmp;
            }
            // normalize pivot row
            Rational pivotVal = aug[col][col];
            for (int c = col; c < m; c++) aug[col][c] = aug[col][c].div(pivotVal);
            // eliminate other rows
            for (int r = 0; r < n; r++) {
                if (r == col) continue;
                Rational factor = aug[r][col];
                if (factor.isZero()) continue;
                for (int c = col; c < m; c++) {
                    aug[r][c] = aug[r][c].sub(factor.mul(aug[col][c]));
                }
            }
        }
        Rational[] sol = new Rational[n];
        for (int i = 0; i < n; i++) sol[i] = aug[i][n];
        return sol;
    }

    // Build augmented matrix (Vandermonde) for given k points and solve for coefficients
    static Rational[] solveVandermonde(List<BigInteger[]> points) {
        int k = points.size();
        Rational[][] aug = new Rational[k][k+1];
        for (int r = 0; r < k; r++) {
            BigInteger x = points.get(r)[0];
            BigInteger pow = BigInteger.ONE;
            for (int c = 0; c < k; c++) {
                aug[r][c] = new Rational(pow);
                pow = pow.multiply(x);
            }
            aug[r][k] = new Rational(points.get(r)[1]); // RHS y
        }
        return solveAugmented(aug);
    }

    // generate all combinations of keysList choose k
    static List<List<Integer>> combinations(List<Integer> keysList, int k) {
        List<List<Integer>> res = new ArrayList<>();
        combineRec(keysList, k, 0, new ArrayList<>(), res);
        return res;
    }
    static void combineRec(List<Integer> arr, int k, int idx, List<Integer> cur, List<List<Integer>> out) {
        if (cur.size() == k) { out.add(new ArrayList<>(cur)); return; }
        for (int i = idx; i < arr.size(); i++) {
            cur.add(arr.get(i));
            combineRec(arr, k, i+1, cur, out);
            cur.remove(cur.size()-1);
        }
    }

    public static void main(String[] args) {
        try {
            String filename = (args.length > 0) ? args[0] : "shares2.json";
            JSONParser parser = new JSONParser();
            JSONObject obj = (JSONObject) parser.parse(new FileReader(filename));

            JSONObject keys = (JSONObject) obj.get("keys");
            int n = Integer.parseInt(keys.get("n").toString());
            int k = Integer.parseInt(keys.get("k").toString());

            // Collect shares into map x -> BigInteger(y)
            Map<Integer, BigInteger> shareMap = new TreeMap<>(); // sorted keys
            if (obj.containsKey("shares")) {
                JSONArray arr = (JSONArray) obj.get("shares");
                for (Object o : arr) {
                    JSONObject s = (JSONObject) o;
                    int x = Integer.parseInt(s.get("x").toString());
                    int base = Integer.parseInt(s.get("base").toString());
                    BigInteger y = new BigInteger(s.get("value").toString(), base);
                    shareMap.put(x, y);
                }
            } else {
                // Hashira format numeric keys: iterate through provided keys in JSON (not rely only 1..n)
                List<String> keyNames = new ArrayList<>();
                for (Object ko : obj.keySet()) {
                    String ks = (String) ko;
                    if (ks.equals("keys")) continue;
                    keyNames.add(ks);
                }
                // sort numeric
                Collections.sort(keyNames, Comparator.comparingInt(Integer::parseInt));
                for (String ks : keyNames) {
                    JSONObject s = (JSONObject) obj.get(ks);
                    int x = Integer.parseInt(ks);
                    int base = Integer.parseInt(s.get("base").toString());
                    BigInteger y = new BigInteger(s.get("value").toString(), base);
                    shareMap.put(x, y);
                }
            }

            // quick sanity
            if (shareMap.size() < k) {
                System.out.println("Not enough shares to reconstruct (have " + shareMap.size() + ", need " + k + ")");
                return;
            }

            List<Integer> keysList = new ArrayList<>(shareMap.keySet());

            // Try all combinations
            List<List<Integer>> allComb = combinations(keysList, k);
            int totalComb = allComb.size();

            Map<Rational, Integer> freq = new HashMap<>();
            Map<Rational, List<Set<Integer>>> sets = new HashMap<>();

            for (List<Integer> comb : allComb) {
                // build points for this combo
                List<BigInteger[]> pts = new ArrayList<>();
                for (int x : comb) pts.add(new BigInteger[]{ BigInteger.valueOf(x), shareMap.get(x) });
                try {
                    Rational[] coeffs = solveVandermonde(pts);
                    if (coeffs == null) continue;
                    Rational secret = coeffs[0]; // a0
                    freq.put(secret, freq.getOrDefault(secret, 0) + 1);
                    sets.computeIfAbsent(secret, z -> new ArrayList<>()).add(new HashSet<>(comb));
                } catch (Exception ex) {
                    // skip singular / failed combos
                }
            }

            if (freq.isEmpty()) {
                System.out.println("Failed to reconstruct secret from any combination.");
                return;
            }

            // choose consensus (max frequency). tie-break by rational's string compare
            Rational bestSecret = null;
            int bestCount = 0;
            for (Map.Entry<Rational,Integer> e : freq.entrySet()) {
                int c = e.getValue();
                if (c > bestCount || (c == bestCount && (bestSecret == null || e.getKey().toString().compareTo(bestSecret.toString()) < 0))) {
                    bestCount = c;
                    bestSecret = e.getKey();
                }
            }

            // Collect keys that appear in any set that yielded the consensus
            Set<Integer> validKeys = new TreeSet<>();
            for (Set<Integer> s : sets.get(bestSecret)) validKeys.addAll(s);

            // Print output like screenshot
            System.out.println("Recovered Secret (a0): " + bestSecret.toString());
            System.out.println("Consensus found in " + bestCount + " combinations out of " + totalComb + ".");
            System.out.println();
            for (int key : validKeys) {
                System.out.println("Key: " + key + "   decimal: " + shareMap.get(key).toString());
            }
            // invalid shares
        Set<Integer> wrongKeys = new TreeSet<>(shareMap.keySet());
        wrongKeys.removeAll(validKeys);
        System.out.println();
        System.out.println("invalid shares:");
        for (int key : wrongKeys) {
            System.out.println("key: " + key + "   decimal: " + shareMap.get(key).toString());
        }

        } catch (Exception e) {
            System.err.println("error: " + e.getMessage());
            e.printStackTrace();
        }
    }
}
