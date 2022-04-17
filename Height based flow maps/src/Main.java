import javax.imageio.ImageIO;
import javax.imageio.stream.FileImageOutputStream;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.List;

public class Main extends JFrame implements MouseWheelListener {


    public static int NR_OF_ROWS;
    public static int NR_OF_COLUMNS;

    public static int source_x;
    public static int source_y;

    public static String TARGET_NAME;

    static JFrame jframe;


    private static double zoomFactor = 500;
    private static double prevZoomFactor = 500;
    private static boolean zoomer;


    public static DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm_ss");
    public static LocalDateTime now = LocalDateTime.now();
    public static String storage_location_name = "";
    public static String currentWorkingPath;


    public static double HEIGHT_FUNCTION_WIDTH;
    public static double HEIGHT_FUNCTION_SCALE;
    public static int NR_OF_ITERATIONS;
    public static String BASE_HEIGHT_TYPE;
    public static double BASE_SCALE;

    public static String DISTANCE_METRIC;

    public static String INPUT_FILE_NAME;

    public static Cell source_cell;

    public static void main(String[] args) throws IOException {

        initialize_parameters();

        storage_location_name = dtf.format(now);
        storage_location_name = storage_location_name.concat("_" + DISTANCE_METRIC + "_" + BASE_HEIGHT_TYPE);

        currentWorkingPath = System.getProperty("user.dir").concat("\\experiments\\");

        File dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\"));
        dir.mkdir();


        double[] widths = {90};
        double[] scales = {10000, 20000, 30000};

        for (int i = 0; i < widths.length; i++) {

            double width = widths[i];
            HEIGHT_FUNCTION_WIDTH = width;

            for (int j = 0; j < scales.length; j++) {
                double scale = scales[j];
                HEIGHT_FUNCTION_SCALE = scale;

                String iteration_location = ("w_" + width + "_s_" + scale);

                dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location + "\\"));

                dir.mkdir();

                Cell[][] grid = initialize_grid(NR_OF_ROWS, NR_OF_COLUMNS);

                ArrayList<String> items = read_input();

                ArrayList<Point> points_list = process_input(items);

                Bounds bounds = obtain_bounds(points_list);

                compute_cell_for_point(bounds, points_list);

                initialize_points_in_grid(grid, points_list);


                if (BASE_HEIGHT_TYPE.equals("EUCLID")) {
                    initialize_grid_height_Euclidean_dist(grid, BASE_SCALE);
                } else if (BASE_HEIGHT_TYPE.equals("default")) {
                    initialize_grid_height(grid, BASE_SCALE);
                } else if (BASE_HEIGHT_TYPE.equals("chebyshev")) {
                    initialize_grid_height_chebyshev_distance(grid, BASE_SCALE);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQRT")) {
                    initialize_grid_height_Euclidean_dist_sqrt(grid, BASE_SCALE);
                }

                compute_flow(grid);

                ArrayList<ArrayList<Cell>> paths = compute_paths(points_list, grid);

                draw(grid, paths, 0, false, iteration_location);

                ArrayList<double[][]> distances_for_paths = null;
                if (DISTANCE_METRIC.equals("BFS")) {
                    distances_for_paths = compute_bfs(grid, paths);
                } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                    distances_for_paths = compute_Dijkstra(grid, paths);
                } else if (DISTANCE_METRIC.equals("ANGULAR")) {
                    distances_for_paths = compute_angular_distance_precise(grid, paths);
                }

                Tuple<Cell[][], ArrayList<ArrayList<Cell>>> result = iterate(grid, points_list, paths, distances_for_paths, NR_OF_ITERATIONS,
                        BASE_HEIGHT_TYPE, BASE_SCALE, true, false, width, scale, iteration_location);

                grid = result.first;

                paths = result.second;

                draw(grid, paths, NR_OF_ITERATIONS, false, iteration_location);

                generate_gif(iteration_location);

                write_output_configuration(iteration_location);

            }
        }
    }

    public static void initialize_parameters() {

        NR_OF_ROWS = 500;
        NR_OF_COLUMNS = 500;
        TARGET_NAME = "A";
        INPUT_FILE_NAME = "./input/1_s_2_t.csv";
        BASE_HEIGHT_TYPE = "EUCLID";
        //BASE_HEIGHT_TYPE = "default";
        //BASE_HEIGHT_TYPE = "chebyshev";
        //BASE_HEIGHT_TYPE = "EUCLID_SQRT";
        DISTANCE_METRIC = "DIJKSTRA";
        //DISTANCE_METRIC = "BFS";
        //DISTANCE_METRIC = "ANGULAR";
        //BASE_SCALE = 0.5;
        NR_OF_ITERATIONS = 3;
        //HEIGHT_FUNCTION_WIDTH = 50; //90;
        //HEIGHT_FUNCTION_SCALE = 200000;//1000000.0;//10000.0; //27000000;//200000;

    }

    public static void write_output_configuration(String iteration_location) throws IOException {
        ////        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/image_" + image_index + ".png"));
        FileWriter fileWriter = new FileWriter(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location) + "\\config.txt");
        fileWriter.write("GRID_SIZE = " + NR_OF_ROWS + " rows x " + NR_OF_COLUMNS + " columns" + "\n");
        fileWriter.write("NR_OF_ITERATIONS = " + NR_OF_ITERATIONS + "\n");
        fileWriter.write("BASE_HEIGHT_TYPE = " + BASE_HEIGHT_TYPE + "\n");
        fileWriter.write("BASE_SCALE = " + BASE_SCALE + "\n");
        fileWriter.write("HEIGHT_FUNCTION_SCALE = " + HEIGHT_FUNCTION_SCALE + "\n");
        fileWriter.write("HEIGHT_FUNCTION_WIDTH = " + HEIGHT_FUNCTION_WIDTH + "\n");
        fileWriter.write("DISTANCE_METRIC = " + DISTANCE_METRIC + "\n");
        fileWriter.close();
    }

    public static void generate_gif(String iteration_location) throws IOException {

//        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/image_" + image_index + ".png"));
        File dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location + "\\"));
        File[] files = dir.listFiles((dir1, name) -> name.endsWith(".png"));


        BufferedImage first = ImageIO.read(new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location + "\\image_0.png")));
        ImageOutputStream output = new FileImageOutputStream(new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location + "\\output_gif.gif")));

        GifSequenceWriter writer = new GifSequenceWriter(output, first.getType(), 250, true);
        writer.writeToSequence(first);


        Arrays.sort(files, new Comparator<File>() {
            @Override
            public int compare(File f1, File f2) {
                String s1 = f1.getName().substring(6, f1.getName().indexOf("."));
                String s2 = f2.getName().substring(6, f2.getName().indexOf("."));
                return Integer.valueOf(s1).compareTo(Integer.valueOf(s2));
            }
        });

        for (File image : files) {
            BufferedImage next = ImageIO.read(image);
            writer.writeToSequence(next);
        }

        writer.close();
        output.close();
    }

    public static Tuple<Cell[][], ArrayList<ArrayList<Cell>>> iterate(
            Cell[][] grid, ArrayList points_list, ArrayList paths,
            ArrayList distances_for_paths, int NR_OF_ITERATIONS,
            String BASE_HEIGHT_TYPE, double base_function_scale, boolean save_outputs,
            boolean show_intermediate_results, double width, double scale, String iteration_location)
            throws IOException {

        adjust_height(grid, distances_for_paths, width, scale);

        // Iterate a number of times:
        for (int i = 1; i < NR_OF_ITERATIONS; i++) {
            System.out.println("iteration : " + i);


            compute_flow(grid);

            paths = compute_paths(points_list, grid);

            if (save_outputs == true) {
                draw(grid, paths, i, show_intermediate_results, iteration_location);
            }

            long startTime = System.currentTimeMillis();

            if (DISTANCE_METRIC.equals("BFS")) {
                distances_for_paths = compute_bfs(grid, paths);
            } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                distances_for_paths = compute_Dijkstra(grid, paths);
            } else if (DISTANCE_METRIC.equals("ANGULAR")) {
                distances_for_paths = compute_angular_distance_precise(grid, paths);
            }
            long endTime = System.currentTimeMillis();
            System.out.println("That took " + (endTime - startTime) + " milliseconds");

            if (BASE_HEIGHT_TYPE.equals("EUCLID")) {
                initialize_grid_height_Euclidean_dist(grid, base_function_scale);
            } else if (BASE_HEIGHT_TYPE.equals("default")) {
                initialize_grid_height(grid, base_function_scale);
            } else if (BASE_HEIGHT_TYPE.equals("chebyshev")) {
                initialize_grid_height_chebyshev_distance(grid, base_function_scale);
            } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQRT")) {
                initialize_grid_height_Euclidean_dist_sqrt(grid, base_function_scale);
            }

            adjust_height(grid, distances_for_paths, width, scale);

        }

        compute_flow(grid);
        paths = compute_paths(points_list, grid);

        Tuple<Cell[][], ArrayList<ArrayList<Cell>>> tuple = new Tuple<Cell[][], ArrayList<ArrayList<Cell>>>(grid, paths);

        return tuple;

    }

    public static double gaussian(double x, double mu, double sigma) {

        return (double) (1.0 / (Math.sqrt(2.0 * Math.PI) * sigma) * Math.exp(-Math.pow((x - mu) / sigma, 2.0) / 2));

    }

    public static void adjust_height(Cell[][] grid, ArrayList distances_for_paths, double width, double scale) {

        System.out.println("adjusting height");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                Cell cell = grid[j][i];

                Iterator path_iterator = distances_for_paths.iterator();

                ArrayList distances_for_cell = new ArrayList();


                while (path_iterator.hasNext()) {

                    double[][] distances = (double[][]) path_iterator.next();

                    distances_for_cell.add(distances[cell.cell_y][cell.cell_x]);

                }

                double sum = 0;
                for (int k = 0; k < distances_for_cell.size(); k++) {
                    sum = sum + gaussian((double) distances_for_cell.get(k), 0, HEIGHT_FUNCTION_WIDTH);
                }

                double height = -HEIGHT_FUNCTION_SCALE * sum;
                grid[j][i].height = grid[j][i].height + height;

            }
        }

    }

    public static void initialize_grid_height(Cell[][] grid, double scale) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (0.05 * (Math.pow(grid[j][i].cell_x - source_x, 2) + Math.pow(grid[j][i].cell_y - source_y, 2)));

            }
        }
    }

    public static void initialize_grid_height_Euclidean_dist_sqrt(Cell[][] grid, double scale) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (0.05 * Math.sqrt(Math.sqrt(Math.pow(grid[j][i].cell_x - source_x, 2) + Math.pow(grid[j][i].cell_y - source_y, 2))));

            }
        }

    }

    public static void initialize_grid_height_chebyshev_distance(Cell[][] grid, double scale) {

        // Direction vectors
        int dRow[] = {-1, 0, 1, 0};
        int dCol[] = {0, 1, 0, -1};

        boolean[][] visited = new boolean[NR_OF_COLUMNS][NR_OF_ROWS];

        double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                distances[i][j] = 0.0f;

            }
        }

        Queue<Cell> queue = new LinkedList<>();


        Cell source_cell = grid[source_y][source_x];

        // add all cells of a path to queue
        queue.add(source_cell);

        visited[source_cell.cell_y][source_cell.cell_x] = true;

        while (!queue.isEmpty()) {

            Cell cell = queue.peek();

            int x = cell.cell_x;
            int y = cell.cell_y;

            queue.remove();

            for (int i = 0; i < 4; i++) {

                int adj_x = x + dCol[i];
                int adj_y = y + dRow[i];

                if (isValid(visited, adj_y, adj_x)) {

                    queue.add(grid[adj_x][adj_y]);
                    visited[adj_y][adj_x] = true;

                    distances[adj_y][adj_x] = distances[y][x] + 1;

                }
            }
        }

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = distances[j][i];

            }
        }
    }


    public static void initialize_grid_height_Euclidean_dist(Cell[][] grid, double scale) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (0.05 * Math.sqrt(Math.pow(grid[j][i].cell_x - source_x, 2) + Math.pow(grid[j][i].cell_y - source_y, 2)));

            }
        }

    }

    public static void assign_zeros_to_grid(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[i][j].height = (double) 0;

            }
        }
    }


    public static class distanceComparator implements Comparator<Cell> {
        public int compare(Cell a, Cell b) {
            if (a.distance < b.distance) {
                return -1;
            } else if (a.distance > b.distance) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    public static boolean is_inside_grid(int x, int y) {

        if (x >= 0 && x < NR_OF_COLUMNS && y >= 0 && y < NR_OF_ROWS) {
            return true;
        } else {
            return false;
        }
    }

    public static ArrayList compute_Dijkstra(Cell[][] grid, ArrayList paths) {

        System.out.println("Computing Dijkstra");

        Iterator path_iterator = paths.iterator();

        double sqrt_2 = Math.sqrt(2);
        double sqrt_5 = Math.sqrt(5);

//        int dRow[] = {-1, -1, 0, 1, 1, 1, 0, -1};
//        int dCol[] = {0, 1, 1, 1, 0, -1, -1, -1};
//        double weights[] = {1, sqrt_2, 1, sqrt_2, 1, sqrt_2, 1, sqrt_2};

        // horse moves:

        int dRow[] = {-1, -1, 0, 1, 1, 1, 0, -1, -2, -1, 1, 2, 2, 1, -1, -2}; // y
        int dCol[] = {0, 1, 1, 1, 0, -1, -1, -1, 1, 2, 2, 1, -1, -2, -2, -1}; // x

        double weights[] = {1, sqrt_2, 1, sqrt_2, 1, sqrt_2, 1, sqrt_2,
                sqrt_5, sqrt_5, sqrt_5, sqrt_5, sqrt_5, sqrt_5, sqrt_5, sqrt_5};

        // for each path we will have a matrix, in which each cell contains distance to the closest point of a path (for n paths)
        ArrayList distances_for_paths = new ArrayList();

        while (path_iterator.hasNext()) {

            ArrayList path = (ArrayList) path_iterator.next();

            // store results for a single path here
            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // Initialize distances with max distance values
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {
                    distances[j][i] = Integer.MAX_VALUE;
                    grid[j][i].distance = Integer.MAX_VALUE;
                }
            }

            PriorityQueue<Cell> Q = new PriorityQueue<Cell>(NR_OF_COLUMNS * NR_OF_ROWS,
                    new distanceComparator());

            // Loop over each cell on a path:
            Iterator cell_iterator = path.iterator();
            while (cell_iterator.hasNext()) {

                // A cell on a path:
                Cell cell = (Cell) cell_iterator.next();

                //cell.distance = 0.0;
                distances[cell.cell_x][cell.cell_y] = 0.0;
                cell.distance = 0.0;

                Q.add(cell);

            }

            while (!Q.isEmpty()) {

                //Tuple<Cell, Double> nd = Q.poll();
                Cell cell = Q.poll();

                int x = cell.cell_x;
                int y = cell.cell_y;

                for (int i = 0; i < 16; i++) {

                    int adj_x = x + dCol[i];
                    int adj_y = y + dRow[i];
                    double weight = weights[i];

                    if (is_inside_grid(adj_x, adj_y)) {

                        if (distances[adj_x][adj_y] > distances[x][y] + weight) {

                            // If Cell is already been reached once,
                            // remove it from priority queue
                            if (distances[adj_x][adj_y] != Integer.MAX_VALUE) {
                                Cell adj = grid[adj_x][adj_y];//new Cell(rows, cols, dist[rows][cols]);
                                adj.distance = distances[adj_x][adj_y];
                                Q.remove(adj);

                            }

                            // Insert cell with updated distance
                            distances[adj_x][adj_y] = distances[cell.cell_x][cell.cell_y] + weight;

                            grid[adj_x][adj_y].distance = distances[adj_x][adj_y];
                            Q.add(grid[adj_x][adj_y]); //new Cell(rows, cols, dist[rows][cols]));

                        }
                    }
                }
            }

            distances_for_paths.add(transposeMatrix(distances));
        }

        return distances_for_paths;
    }

    public static double[][] transposeMatrix(double[][] matrix) {
        // TODO: figure out how to assign these without transpose
        int m = matrix.length;
        int n = matrix[0].length;

        double[][] transposedMatrix = new double[n][m];

        for (int x = 0; x < n; x++) {
            for (int y = 0; y < m; y++) {
                transposedMatrix[x][y] = matrix[y][x];
            }
        }

        return transposedMatrix;
    }

    public static ArrayList compute_angular_distance_precise(Cell[][] grid, ArrayList paths) {

        System.out.println("computing angular distance");

        Iterator path_iterator = paths.iterator();

        ArrayList distances_for_paths = new ArrayList();

        // for all paths
        while (path_iterator.hasNext()) {

            ArrayList path = (ArrayList) path_iterator.next();

            Collections.reverse(path);

            // the first cell in the path is the source node (A). The last cell is the target (B)
            ArrayList<IntermediateCell> finer_path = new ArrayList<IntermediateCell>();

            //double distance = 0.0;
            // Make path finer:

            double sum = 0.0;

            for (int i = 0; i < path.size() - 1; i++) { // potentially up to size - 1

                Cell cell = (Cell) path.get(i);
                Cell next_cell = (Cell) path.get(i + 1);

//                double distance = (Math.sqrt(Math.pow(next_cell.cell_x - cell.cell_x, 2) +
//                        Math.pow(next_cell.cell_y - cell.cell_y, 2)));

                int cell_x = cell.cell_x;
                int next_cell_x = next_cell.cell_x;

                int cell_y = cell.cell_y;
                int next_cell_y = next_cell.cell_y;

//                double distance_x = Math.sqrt(Math.pow(next_cell.cell_x - cell.cell_x, 2));
//                double distance_y = Math.sqrt(Math.pow(next_cell.cell_y - cell.cell_y, 2));

                double abs_x = Math.abs(cell_x - next_cell_x);
                double abs_y = Math.abs(cell_y - next_cell_y);

                double splits = 100;

                double split_x = abs_x / splits;
                double split_y = abs_y / splits;

                for (int j = 0; j < splits; j++) {

                    IntermediateCell intermediate_cell = new IntermediateCell();

                    intermediate_cell.cell_x = cell_x - split_x * j;
                    intermediate_cell.cell_y = cell_y - split_y * j;

                    finer_path.add(intermediate_cell);
                    //sum = sum + split_distance;
                    //finer_path.add(sum);
                }
            }
            IntermediateCell last_cell = new IntermediateCell();

            Cell last_path_cell = (Cell) path.get(path.size() - 1);
            last_cell.cell_x = last_path_cell.cell_x;
            last_cell.cell_y = last_path_cell.cell_y;

            finer_path.add(last_cell);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_x - cell.cell_x, 2) +
                            Math.pow(source_cell.cell_y - cell.cell_y, 2)));

                    // We don't actually need an index of a cell. We just need the distance
                    // For our binary search: we split the actual path into a much finer path (each cell split into 10)
                    // Then run binary search on the finer path.
                    // Input: Arraylist/array of bins. |array| = number of cells + cumber of finer cells
                    // Each smaller cell holds the distance it corresponds to in the split

                    //double distance = binary_serach_for_finer_path(finer_path, radius);

                    int index_of_cell = binary_serach_for_finer_path(finer_path, radius);

                    double x = finer_path.get(index_of_cell).cell_x;
                    double y = finer_path.get(index_of_cell).cell_y;

                    IntermediateCell intersection_cell = (IntermediateCell) finer_path.get(index_of_cell);

                    // Here we can either use radius as dist or the actual distance. Matter of precision.
                    double dist = (Math.sqrt(Math.pow(intersection_cell.cell_x - source_cell.cell_x, 2) +
                            Math.pow(intersection_cell.cell_y - source_cell.cell_y, 2)));

                    double dist_2 = (Math.sqrt(Math.pow(intersection_cell.cell_x - cell.cell_x, 2) +
                            Math.pow(intersection_cell.cell_y - cell.cell_y, 2)));

                    double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(dist_2, 2)) /
                            (2.0 * radius * dist));

                    distances[i][j] = angle;
                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_x][cell.cell_y] = 0.0;

            }

            distances_for_paths.add(transposeMatrix(distances));
        }
        return distances_for_paths;

    }

    public static ArrayList compute_angular_distance(Cell[][] grid, ArrayList paths) {

        System.out.println("computing angular distance");

        Iterator path_iterator = paths.iterator();

        ArrayList distances_for_paths = new ArrayList();

        // for all paths
        while (path_iterator.hasNext()) {

            ArrayList path = (ArrayList) path_iterator.next();

            Collections.reverse(path);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_x - cell.cell_x, 2) +
                            Math.pow(source_cell.cell_y - cell.cell_y, 2)));

                    int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                    Cell intersection_cell = (Cell) path.get(index_of_cell);

                    // Here we can either use radius as dist or the actual distance. Matter of precision.
                    double dist = (Math.sqrt(Math.pow(intersection_cell.cell_x - source_cell.cell_x, 2) +
                            Math.pow(intersection_cell.cell_y - source_cell.cell_y, 2)));

                    double dist_2 = (Math.sqrt(Math.pow(intersection_cell.cell_x - cell.cell_x, 2) +
                            Math.pow(intersection_cell.cell_y - cell.cell_y, 2)));

                    double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(dist_2, 2)) /
                            (2.0 * radius * dist));

                    distances[i][j] = angle;
                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_x][cell.cell_y] = 0.0;

            }

            distances_for_paths.add(transposeMatrix(distances));
        }
        return distances_for_paths;

    }

    public static int binary_serach_for_finer_path(ArrayList path, double distance_target) {

        int start = 0;
        int end = path.size() - 1;

        int ans = -1;

        while (start <= end) {
            int mid = (start + end) / 2;
            IntermediateCell mid_cell = (IntermediateCell) path.get(mid);

            double dist_to_mid = (Math.sqrt(Math.pow(source_cell.cell_x - mid_cell.cell_x, 2) +
                    Math.pow(source_cell.cell_y - mid_cell.cell_y, 2)));

            if (dist_to_mid == distance_target) {
                return mid;
            }

            if (dist_to_mid <= distance_target) {
                start = mid + 1;
            } else {
                ans = mid;
                end = mid - 1;
            }
        }
        return end;
    }

    public static int binary_search_2(ArrayList<Cell> path, double distance_target) {

        int start = 0;
        int end = path.size() - 1;

        int ans = -1;

        while (start <= end) {
            int mid = (start + end) / 2;
            Cell mid_cell = path.get(mid);

            double dist_to_mid = (Math.sqrt(Math.pow(source_cell.cell_x - mid_cell.cell_x, 2) +
                    Math.pow(source_cell.cell_y - mid_cell.cell_y, 2)));

            if (dist_to_mid == distance_target) {
                return mid;
            }

            if (dist_to_mid <= distance_target) {
                start = mid + 1;
            } else {
                ans = mid;
                end = mid - 1;
            }
        }
        return end;
    }

    public static int binarySearch(ArrayList<Cell> path, int l, int r, double dist) {

        if (r >= l) {
            int mid = l + (r - l) / 2;

            Cell mid_cell = path.get(mid);

            double dist_to_mid = (Math.sqrt(Math.pow(source_cell.cell_x - mid_cell.cell_x, 2) +
                    Math.pow(source_cell.cell_y - mid_cell.cell_y, 2)));

            if (dist_to_mid == dist) {
                return mid;
            } else if (dist_to_mid > dist) {
                binarySearch(path, l, mid - 1, dist);
            } else {
                binarySearch(path, mid + 1, r, dist);
            }
        }

        return -1;
    }

    public static ArrayList compute_angular_distance_2(Cell[][] grid, ArrayList paths) {

        System.out.println("computing angular distance");

        Iterator path_iterator = paths.iterator();

        ArrayList distances_for_paths = new ArrayList();

        Cell source_cell = grid[source_x][source_y];

        // for all paths
        while (path_iterator.hasNext()) {

            ArrayList path = (ArrayList) path_iterator.next();

            Collections.reverse(path);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    //System.out.println("i:" + i + " j: " + j);
                    Cell cell = grid[i][j];

                    if (path.contains(cell)) {
                        distances[i][j] = 0.0;
                        continue;
                    }

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_x - cell.cell_x, 2) +
                            Math.pow(source_cell.cell_y - cell.cell_y, 2)));

                    Iterator cell_iterator = path.iterator();

                    while (cell_iterator.hasNext()) {

                        Cell path_cell = (Cell) cell_iterator.next();

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_x - path_cell.cell_x, 2) +
                                Math.pow(source_cell.cell_y - path_cell.cell_y, 2)));
                        if (dist <= radius) {
                            // consider next cell on path
                            continue;
                        } else {
                            // this is the cell
                            double dist_2 = (Math.sqrt(Math.pow(cell.cell_x - path_cell.cell_x, 2) +
                                    Math.pow(cell.cell_y - path_cell.cell_y, 2)));

                            double angle = Math.acos((Math.pow(dist, 2) + Math.pow(radius, 2) - Math.pow(dist_2, 2)) /
                                    (2.0 * dist * radius));

                            distances[i][j] = angle;

                        }
                    }
                }
            }

            distances_for_paths.add(transposeMatrix(distances));
        }
        return distances_for_paths;
    }

    public static ArrayList compute_bfs(Cell[][] grid, ArrayList paths) {

        System.out.println("computing bfs");

        Iterator path_iterator = paths.iterator();


        // Direction vectors
        int dRow[] = {-1, 0, 1, 0};
        int dCol[] = {0, 1, 0, -1};

        // dw = weight of the edge

        ArrayList distances_for_paths = new ArrayList();

        while (path_iterator.hasNext()) {


            boolean[][] visited = new boolean[NR_OF_COLUMNS][NR_OF_ROWS];

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    distances[i][j] = 0.0f;

                }
            }

            ArrayList path = (ArrayList) path_iterator.next();

            Queue<Cell> queue = new LinkedList<>();

            Iterator cell_iterator = path.iterator();

            while (cell_iterator.hasNext()) {

                Cell cell = (Cell) cell_iterator.next();

                // add all cells of a path to queue
                queue.add(cell);

                visited[cell.cell_y][cell.cell_x] = true;

            }

            while (!queue.isEmpty()) {

                Cell cell = queue.peek();

                int x = cell.cell_x;
                int y = cell.cell_y;

                queue.remove();

                for (int i = 0; i < 4; i++) {

                    int adj_x = x + dCol[i];
                    int adj_y = y + dRow[i];

                    if (isValid(visited, adj_y, adj_x)) {

                        queue.add(grid[adj_x][adj_y]);
                        visited[adj_y][adj_x] = true;

                        distances[adj_y][adj_x] = distances[y][x] + 1;
                    }
                }
            }

            distances_for_paths.add(distances);
        }

        return distances_for_paths;
    }


    static boolean isValid(boolean visited[][], int row, int col) {

        // If cell lies out of bounds
        if (row < 0 || col < 0 || row >= NR_OF_ROWS || col >= NR_OF_COLUMNS)
            return false;

        // If cell is already visited
        if (visited[row][col])
            return false;

        // Otherwise
        return true;
    }

    public static void compute_shortest_paths_naive(Cell[][] grid, ArrayList paths) {

        //https://programmer.ink/think/graph-theory-search-how-to-use-multi-source-bfs-to-reduce-time-complexity.html
        //https://www.geeksforgeeks.org/multi-source-shortest-path-in-unweighted-graph/

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                Cell current_cell = grid[i][j];

                Iterator path_it = paths.iterator();

                while (path_it.hasNext()) {

                    ArrayList path = (ArrayList) path_it.next();

                    Iterator path_iter = path.iterator();

                    ArrayList distances = new ArrayList();

                    while (path_iter.hasNext()) {

                        Cell path_cell = (Cell) path_iter.next();

                        double distance = (double) (Math.sqrt(Math.pow(path_cell.cell_x - current_cell.cell_x, 2) + Math.pow(path_cell.cell_y - current_cell.cell_y, 2)));

                        distances.add(distance);

                    }

                    int min_distance_index = distances.indexOf(Collections.min(distances));

                    Cell min_distance_cell = (Cell) path.get(min_distance_index);

                }
            }
        }
    }

    public static void draw_distances(Cell[][] grid, ArrayList paths, ArrayList distances_for_paths) throws
            IOException {

        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        double[][] dist = (double[][]) distances_for_paths.get(1);

        double max_dist = dist[0][0];
        double min_dist = 0;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                if (dist[i][j] > max_dist) {
                    max_dist = dist[i][j];
                }
            }
        }

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                int c = (int) (dist[i][j] * 255.0 / max_dist);


                if (c < 0) {
                    image.setRGB(i, j, new Color(0, 0, -c).getRGB());
                } else {
                    image.setRGB(i, j, new Color(c, 0, 0).getRGB());
                }


//                if (!grid[i][j].title.equals("")) {
//
//                    image.setRGB(i, j, new Color(0, 255, 0).getRGB());
//                }
            }
        }

        Iterator iter = paths.iterator();

        while (iter.hasNext()) {

            ArrayList path = (ArrayList) iter.next();

            Iterator cell_iter = path.iterator();

            while (cell_iter.hasNext()) {

                Cell cell = (Cell) cell_iter.next();

                image.setRGB((int) cell.cell_x, (int) cell.cell_y, new Color(255, 255, 255).getRGB());

            }
        }

        JPanel pane = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                g.drawImage(image, 0, 0, null);
            }
        };

        jframe.add(pane);

        jframe.setVisible(true);
        jframe.show();

        ImageIO.write(image, "png", new File("image.png"));

    }

    public static Color getValueBetweenTwoFixedColors(float value) {
        int aR = 0;
        int aG = 0;
        int aB = 255;  // RGB for our 1st color (blue in this case).
        int bR = 255;
        int bG = 0;
        int bB = 0;    // RGB for our 2nd color (red in this case).

        int red = (int) ((float) (bR - aR) * value + aR);      // Evaluated as -255*value + 255.
        int green = (int) ((float) (bG - aG) * value + aG);      // Evaluates as 0.
        int blue = (int) ((float) (bB - aB) * value + aB);      // Evaluates as 255*value + 0.

        Color color = new Color(red, green, blue);

        return color;
    }

    public static Color getHeatMapColor(float value) {

        // TODO: fix, low priority
        int NUM_COLORS = 4;

        float color[][] = {{0, 0, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}};
        // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.

        int idx1;        // |-- Our desired color will be between these two indexes in "color".
        int idx2;        // |
        float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

        if (value <= 0) {
            idx1 = idx2 = 0;
        }    // accounts for an input <=0
        else if (value >= 1) {
            idx1 = idx2 = NUM_COLORS - 1;
        }    // accounts for an input >=0
        else {
            value = value * (NUM_COLORS - 1);        // Will multiply value by 3.
            idx1 = (int) Math.floor(value);                  // Our desired color will be after this index.
            idx2 = idx1 + 1;                        // ... and before this index (inclusive).
            fractBetween = value - idx1;    // Distance between the two indexes (0-1).
        }

        int red = (int) ((color[idx2][0] - color[idx1][0]) * fractBetween + color[idx1][0]);
        int green = (int) ((color[idx2][1] - color[idx1][1]) * fractBetween + color[idx1][1]);
        int blue = (int) ((color[idx2][2] - color[idx1][2]) * fractBetween + color[idx1][2]);

        Color result_color = new Color(255 * red, 255 * green, 255 * blue);

        return result_color;
    }

    public static void draw(Cell[][] grid, ArrayList paths, int image_index, boolean show_intermediate_results, String iteration_location)
            throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        double max_height = grid[0][0].height;
        double min_hieght = grid[0][0].height;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                if (grid[i][j].height > max_height) {
                    max_height = grid[i][j].height;
                }
                if (grid[i][j].height < min_hieght) {
                    min_hieght = grid[i][j].height;
                }
            }
        }

        System.out.println("min : " + min_hieght + " max: " + max_height);

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                float value = (float) ((grid[i][j].height - min_hieght) / (max_height - min_hieght));

                Color color = getValueBetweenTwoFixedColors(value);
//                Color colorOfYourDataPoint = null;
//                try {
//                    colorOfYourDataPoint = new Color((int) value*255, (int)value*255, (int)value*255);
//                } catch (Exception e ) {
//                    System.out.println();
//                }

                image.setRGB(i, j, color.getRGB());

                //int c = (int) (grid[i][j].height * 255.0 / max_height);
//                if (c < 0) {
//                    c = (int) (grid[i][j].height * 255.0 / min_hieght);
//                    image.setRGB(i, j, new Color(0, 0, c).getRGB());
//                } else {
//                    c = (int) (grid[i][j].height * 255.0 / max_height);
//                    try {
//                        image.setRGB(i, j, new Color(c, 0, 0).getRGB());
//                    } catch (Exception e) {
//                        System.out.println();
//                    }
//                }


                if (!grid[i][j].title.equals("")) {

                    image.setRGB(i, j, new Color(0, 255, 0).getRGB());
                }
            }
        }

        Iterator iter = paths.iterator();

        while (iter.hasNext()) {

            ArrayList path = (ArrayList) iter.next();

            Iterator cell_iter = path.iterator();

            while (cell_iter.hasNext()) {

                Cell cell = (Cell) cell_iter.next();

                image.setRGB((int) cell.cell_x, (int) cell.cell_y, new Color(255, 255, 255).getRGB());

            }
        }


        if (show_intermediate_results) {

            JPanel pane = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    g.drawImage(image, 0, 0, null);
                }
            };

            jframe.add(pane);

            jframe.setVisible(true);

            jframe.show();
        }
        //dir = new File(currentWorkingPath.concat("\\" +  storage_location_name + "\\" + iteration_location + "\\"));

        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/image_" + image_index + ".png"));
        file.mkdirs();
        ImageIO.write(image, "png", file);

    }

    @Override
    public void paint(Graphics g) {
        super.paint(g);
        Graphics2D g2 = (Graphics2D) g;
        if (zoomer) {
            AffineTransform at = new AffineTransform();
            at.scale(zoomFactor, zoomFactor);
            prevZoomFactor = zoomFactor;
            g2.transform(at);
            zoomer = false;
        }
        // All drawings go here
    }


    public static ArrayList compute_paths(ArrayList points_list, Cell[][] grid) {

        System.out.println("computing paths");

        Iterator it = points_list.iterator();

        ArrayList<ArrayList<Cell>> paths = new ArrayList();

        while (it.hasNext()) {

            Point point = (Point) it.next();

            if (point.name.equals(TARGET_NAME)) {
                continue;
            }

            Cell grid_cell = grid[point.grid_x][point.grid_y];

            Cell current_cell = grid_cell;

            ArrayList<Cell> path = new ArrayList();

            path.add(current_cell);

            while (!(current_cell.title.equals(TARGET_NAME))) {

                int node_x = (int) current_cell.cell_x;
                int node_y = (int) current_cell.cell_y;

                double flow = current_cell.flow_direction;

                if (flow == 1) {
                    current_cell = grid[node_x + 1][node_y];
                } else if (flow == 2) {
                    current_cell = grid[node_x + 1][node_y + 1];
                } else if (flow == 4) {
                    current_cell = grid[node_x][node_y + 1];
                } else if (flow == 8) {
                    current_cell = grid[node_x - 1][node_y + 1];
                } else if (flow == 16) {
                    current_cell = grid[node_x - 1][node_y];
                } else if (flow == 32) {
                    current_cell = grid[node_x - 1][node_y - 1];
                } else if (flow == 64) {
                    current_cell = grid[node_x][node_y - 1];
                } else if (flow == 128) {
                    current_cell = grid[node_x + 1][node_y - 1];
                }

                path.add(current_cell);

            }

            paths.add(path);

        }
        return paths;
    }

    public static void compute_flow(Cell[][] grid) {
        System.out.println("computing flow");
        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                int x = (int) grid[i][j].cell_x;
                int y = (int) grid[i][j].cell_y;

                ArrayList<Cell> neighbors = new ArrayList();

                Cell left = null;
                Cell right = null;
                Cell top = null;
                Cell bottom = null;
                Cell top_left = null;
                Cell top_right = null;
                Cell bottom_left = null;
                Cell bottom_right = null;

                if (x - 1 >= 0) {
                    left = grid[x - 1][y];
                    neighbors.add(left);
                }

                if (y + 1 < NR_OF_ROWS) {
                    bottom = grid[x][y + 1];
                    neighbors.add(bottom);
                }

                if (x + 1 < NR_OF_COLUMNS) {
                    right = grid[x + 1][y];
                    neighbors.add(right);
                }

                if (y - 1 >= 0) {
                    top = grid[x][y - 1];
                    neighbors.add(top);
                }

                if (x - 1 >= 0 && y - 1 >= 0) {
                    top_left = grid[x - 1][y - 1];
                    neighbors.add(top_left);
                }

                if (y - 1 >= 0 && x + 1 < NR_OF_COLUMNS) {
                    top_right = grid[x + 1][y - 1];
                    neighbors.add(top_right);
                }

                if (x - 1 >= 0 && y + 1 < NR_OF_ROWS) {
                    bottom_left = grid[x - 1][y + 1];
                    neighbors.add(bottom_left);
                }

                if (x + 1 < NR_OF_COLUMNS && y + 1 < NR_OF_ROWS) {
                    bottom_right = grid[x + 1][y + 1];
                    neighbors.add(bottom_right);
                }

                Iterator it = neighbors.iterator();

                ArrayList drop_for_neighbors = new ArrayList();

                while (it.hasNext()) {

                    Cell neighbor = (Cell) it.next();

                    double change_in_height = grid[x][y].height - neighbor.height;

                    double distance = 0f;

                    if (neighbor == left || neighbor == right || neighbor == top || neighbor == bottom) {

                        distance = 1f;

                    } else if (neighbor == top_left || neighbor == top_right || neighbor == bottom_left || neighbor == bottom_right) {
                        distance = (double) Math.sqrt(2);
                    }

                    double drop = (change_in_height / distance);
                    drop_for_neighbors.add(drop);

                }

                int max_drop_index = drop_for_neighbors.indexOf(Collections.max(drop_for_neighbors));

                Cell max_drop_neighbor = (Cell) neighbors.get(max_drop_index);

                if (max_drop_neighbor == left) {
                    grid[x][y].flow_direction = 16;
                } else if (max_drop_neighbor == right) {
                    grid[x][y].flow_direction = 1;
                } else if (max_drop_neighbor == bottom) {
                    grid[x][y].flow_direction = 4;
                } else if (max_drop_neighbor == top) {
                    grid[x][y].flow_direction = 64;
                } else if (max_drop_neighbor == top_left) {
                    grid[x][y].flow_direction = 32;
                } else if (max_drop_neighbor == top_right) {
                    grid[x][y].flow_direction = 128;
                } else if (max_drop_neighbor == bottom_left) {
                    grid[x][y].flow_direction = 8;
                } else if (max_drop_neighbor == bottom_right) {
                    grid[x][y].flow_direction = 2;
                }

            }
        }
    }


    public static void initialize_points_in_grid(Cell[][] grid, ArrayList<Point> points_list) {

        Iterator it = points_list.iterator();

        // setting the names of source/targets and setting height of source
        while (it.hasNext()) {

            Point point = (Point) it.next();

            grid[point.grid_x][point.grid_y].title = point.name;

            if (point.name.equals(TARGET_NAME)) {
                grid[point.grid_x][point.grid_y].height = 0;

                source_x = point.grid_x;
                source_y = point.grid_y;

            }
        }

        source_cell = grid[source_x][source_y];
    }


    public static void compute_cell_for_point(Bounds bounds, ArrayList points) {

        Iterator it = points.iterator();

        double extent_x = bounds.max_x - bounds.min_x;
        double extent_y = bounds.max_y - bounds.min_y;

        double step_x = extent_x / NR_OF_COLUMNS;
        double step_y = extent_y / NR_OF_ROWS;

        while (it.hasNext()) {

            Point point = (Point) it.next();

            int i = 1;
            int j = 1;

            while (bounds.min_x + i * step_x < point.x) {
                i = i + 1;
            }
            int column_index = i - 1;

            while (bounds.min_y + j * step_y < point.y) {
                j = j + 1;
            }
            int row_index = j - 1;

            point.grid_x = column_index;
            point.grid_y = row_index;

        }

    }

    public static ArrayList process_input(ArrayList items) {

        Iterator it = items.iterator();

        ArrayList point_list = new ArrayList();

        while (it.hasNext()) {

            String input_line = it.next().toString();
            List<String> components = Arrays.asList(input_line.split(";"));

            Point point = new Point();
            point.name = components.get(0);
            point.x = Float.parseFloat(components.get(1));
            point.y = Float.parseFloat(components.get(2));

            point_list.add(point);

        }

        return point_list;
    }

    public static Bounds obtain_bounds(ArrayList<Point> input_points) {

        Iterator it = input_points.iterator();

        double min_x = input_points.get(0).x;
        double min_y = input_points.get(0).y;
        double max_x = input_points.get(0).x;
        double max_y = input_points.get(0).y;

        while (it.hasNext()) {

            Point next_point = (Point) it.next();

            if (next_point.x < min_x) {
                min_x = next_point.x;
            }
            if (next_point.x > max_x) {
                max_x = next_point.x;
            }
            if (next_point.y < min_y) {
                min_y = next_point.y;
            }
            if (next_point.y > max_y) {
                max_y = next_point.y;
            }

        }

        Bounds result = new Bounds();
        result.max_x = max_x;
        result.min_x = min_x;
        result.min_y = min_y;
        result.max_y = max_y;

        return result;
    }


    public static ArrayList read_input() throws FileNotFoundException {

        Scanner sc = new Scanner(new File(INPUT_FILE_NAME));
        sc.useDelimiter("\n");

        ArrayList items = new ArrayList();

        while (sc.hasNext()) {

            items.add(sc.next());

        }
        sc.close();

        return items;
    }

    public static Cell[][] initialize_grid(int nr_of_rows, int nr_of_columns) {

        Cell[][] grid = new Cell[nr_of_rows][nr_of_columns];

        for (int i = 0; i < nr_of_rows; i++) {
            for (int j = 0; j < nr_of_columns; j++) {

                Cell cell = new Cell();

                cell.cell_x = i;
                cell.cell_y = j;
                cell.flow_direction = 0;
                cell.height = 0;

                cell.title = "";

                grid[i][j] = cell;

            }
        }

        return grid;
    }

    @Override
    public void mouseWheelMoved(MouseWheelEvent mouseWheelEvent) {
        zoomer = true;
        //Zoom in
        if (mouseWheelEvent.getWheelRotation() < 0) {
            zoomFactor *= 1.1;
            repaint();
        }
        //Zoom out
        if (mouseWheelEvent.getWheelRotation() > 0) {
            zoomFactor /= 1.1;
            repaint();
        }
    }
}
