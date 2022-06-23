public class BicubicInterpolator extends CubicInterpolator {
    private double[] arr = new double[4];

    // p - matrix 4x4 in which the middle  2x2 grid is where the interpolation happens
    // x, y - values in grid that we want to interpolate
    public double getValue(double[][] p, double x, double y) {
        arr[0] = getValue(p[0], y);
        arr[1] = getValue(p[1], y);
        arr[2] = getValue(p[2], y);
        arr[3] = getValue(p[3], y);
        return getValue(arr, x);
    }
}
