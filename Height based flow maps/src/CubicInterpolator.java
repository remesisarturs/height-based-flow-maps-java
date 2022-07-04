public class CubicInterpolator {
    public static double getValue(double[] p, double x) {

        double solution = p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0 * p[0] - 5.0 * p[1] + 4.0 * p[2] - p[3] + x * (3.0 * (p[1] - p[2]) + p[3] - p[0])));

        // p_1 + 0.5 * x * (p_2 - p_0 + x * (2.0 * p_0 - 5.0 * p_1 + 4.0 * p_2 - p_3 + x * (3.0 * (p_1 - p_2) + p_3 - p_0)))

        double a = ((-1/2) * p[0]) + ((3/2) * p[1]) - ((3/2) * p[2]) + ((1/2) * p[3]);
        double b = p[0] - ((5/2) * p[1]) + (2 * p[2]) - ((1/2) * p[3]);
        double c = ((-1/2) * p[0]) + ((1/2) * p[2]);
        double d = p[1];

        //double solution_2 = a * Math.pow(x, 3) + b * Math.pow(x, 2) + c * x + d;

        System.out.println();

        return solution;

    }
}
