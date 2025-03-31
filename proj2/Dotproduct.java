import java.util.*;
import java.util.regex.Pattern;

public class Dotproduct {
    enum DataType { CHAR, INT, DOUBLE }

    static class Element {
        final DataType type;
        final Object value;
        
        Element(DataType type, Object value) {
            this.type = type;
            this.value = value;
        }
        
        static Element parse(String s) throws IllegalArgumentException {
            try {
                return new Element(DataType.INT, Integer.parseInt(s));
            } catch (NumberFormatException ignored) {}
            
            try {
                return new Element(DataType.DOUBLE, Double.parseDouble(s));
            } catch (NumberFormatException ignored) {}
            
            if (s.length() == 1) {
                return new Element(DataType.CHAR, s.charAt(0));
            }
            
            throw new IllegalArgumentException("Invalid element: " + s);
        }
    }

    public static void main(String[] args) {
        try {
            Scanner scanner = new Scanner(System.in);
            scanner.useDelimiter("\\R"); 
            
            List<Element> vec1 = readVector(scanner);
            List<Element> vec2 = readVector(scanner);
            scanner.close();

            if (vec1.size() != vec2.size()) {
                throw new IllegalArgumentException(String.format(
                    "Vectors length mismatch (%d vs %d)", vec1.size(), vec2.size()));
            }

            long startTime = System.nanoTime();
            double dotProduct = calculateDotProduct(vec1, vec2);
            double elapsed = (System.nanoTime() - startTime) / 1e9;

            printResult(dotProduct);
            System.out.printf("%nTime elapsed: %.9f seconds%n", elapsed);
            System.out.flush();
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            System.exit(1);
        }
    }

    private static List<Element> readVector(Scanner scanner) {
        String input = scanner.hasNext() ? scanner.next().trim() : "";
        if (input.isEmpty()) {
            throw new IllegalArgumentException("Empty input");
        }

        String[] elements = input.split("[ ,;]+");
        List<Element> vector = new ArrayList<>();
        for (String elem : elements) {
            vector.add(Element.parse(elem));
        }
        return vector;
    }

    private static double calculateDotProduct(List<Element> vec1, List<Element> vec2) {
        double total = 0.0;
        for (int i = 0; i < vec1.size(); i++) {
            Element e1 = vec1.get(i);
            Element e2 = vec2.get(i);
            
            DataType targetType = higherType(e1.type, e2.type);
            double v1 = convertValue(e1, targetType);
            double v2 = convertValue(e2, targetType);
            
            total += v1 * v2;
        }
        return total;
    }

    private static DataType higherType(DataType t1, DataType t2) {
        return t1.ordinal() > t2.ordinal() ? t1 : t2;
    }

    private static double convertValue(Element e, DataType target) {
        switch (e.type) {
            case CHAR:   return convertChar((Character)e.value, target);
            case INT:    return convertInt((Integer)e.value, target);
            default:     return (Double)e.value;
        }
    }

    private static double convertChar(char c, DataType target) {
        return (target == DataType.CHAR) ? c : (target == DataType.INT) ? (int)c : (double)c;
    }

    private static double convertInt(int i, DataType target) {
        return (target == DataType.INT) ? i : (double)i;
    }

    private static void printResult(double result) {
        if (result == (int)result) {
            System.out.printf("Dot product: %d%n", (int)result);
        } else {
            System.out.printf("Dot product: %.2f%n", result);
        }
    }
}