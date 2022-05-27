import javafx.util.Pair;

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


    private static double zoomFactor = 50;
    private static double prevZoomFactor = 50;
    private static boolean zoomer;


    public static DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm_ss");
    public static LocalDateTime now = LocalDateTime.now();
    public static String storage_location_name = "";
    public static String currentWorkingPath;

    public static int GIF_DELAY;


    public static double HEIGHT_FUNCTION_WIDTH;
    public static double HEIGHT_FUNCTION_SCALE;
    public static int NR_OF_ITERATIONS;
    public static String BASE_HEIGHT_TYPE;
    public static double BASE_SCALE;

    public static double ARC_RADIUS;

    public static String DISTANCE_METRIC;

    public static String INPUT_FILE_NAME;

    public static Cell source_cell;

    public static double[] WIDTHS;
    public static double[] SCALES;

    public static boolean REMOVE_DIAGONAL_BIAS = false;
    public static boolean RESET_HEIGHTS = true;

    public static final String ANSI_RESET = "\u001B[0m";
    public static final String ANSI_BLACK = "\u001B[30m";
    public static final String ANSI_RED = "\u001B[31m";
    public static final String ANSI_GREEN = "\u001B[32m";
    public static final String ANSI_YELLOW = "\u001B[33m";
    public static final String ANSI_BLUE = "\u001B[34m";
    public static final String ANSI_PURPLE = "\u001B[35m";
    public static final String ANSI_CYAN = "\u001B[36m";
    public static final String ANSI_WHITE = "\u001B[37m";
    public static boolean DRAW_TEXT_DESCRIPTION = true;
    public static boolean DRAW_PATHS = true;
    public static double previous_sum = 0.0;
    public static double new_sum = 0.0;
    public static boolean DRAW_DISTANCE_IMAGES = false;
    public static boolean HEIGHT_DECAY_ENABLED = false;
    public static double MIN_HEIGHT = Integer.MAX_VALUE;
    public static double MAX_HEIGHT = Integer.MIN_VALUE;
    public static boolean GENERATE_INTERMEDIATE_RESULTS = true;
    public static boolean GENERATE_INTERMEDIATE_HEIGHT = true;
    public static boolean EXPERIMENTAL_MODE = false;
    public static FileWriter log_file_writer;
    public static String COLOR_MODE = "";
    public static boolean PATH_SCALING = false;
    public static boolean MEMORY_MODE = false;

    public static double MEMORY_DECAY_RATE;

    public static String SCALING_MODE;

    public static void main(String[] args) throws IOException {

        initialize_parameters();

        storage_location_name = dtf.format(now);
        storage_location_name = storage_location_name.concat("_" + DISTANCE_METRIC + "_" + BASE_HEIGHT_TYPE);

        currentWorkingPath = System.getProperty("user.dir").concat("\\experiments\\");

        File dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\"));
        dir.mkdir();

        for (int i = 0; i < WIDTHS.length; i++) {

            double width = WIDTHS[i];
            HEIGHT_FUNCTION_WIDTH = width;

            for (int j = 0; j < SCALES.length; j++) {
                double scale = SCALES[j];
                HEIGHT_FUNCTION_SCALE = scale;

                String iteration_location = ("w_" + width + "_s_" + scale);

                dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location + "\\"));

                dir.mkdir();

                log_file_writer = new FileWriter(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location) + "\\log.txt");

                write_output_configuration(iteration_location);

                Cell[][] grid = initialize_grid(NR_OF_ROWS, NR_OF_COLUMNS);
                double[][] base_function = new double[NR_OF_ROWS][NR_OF_COLUMNS];

                ArrayList<String> items = read_input();

                ArrayList<Point> points_list = process_input(items);

                Bounds bounds = obtain_bounds(points_list);

                compute_cell_for_point(bounds, points_list);

                initialize_points_in_grid(grid, points_list);

                if (BASE_HEIGHT_TYPE.equals("EUCLID")) {
                    initialize_grid_height_Euclidean_dist(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQUARED")) {
                    initialize_grid_height_Euclidean_squared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("chebyshev")) {
                    initialize_grid_height_chebyshev_distance(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQRT")) {
                    initialize_grid_height_Euclidean_dist_sqrt(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE")) {
                    initialize_grid_height_to_edge(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQUARED")) {
                    initialize_grid_height_to_edge_squared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQRT")) {
                    initialize_grid_height_to_edge_sqrt(grid);
                }

                compute_flow(grid, 0);

                ArrayList<Path> paths;

                if (EXPERIMENTAL_MODE) {
                    paths = compute_paths_to_frame_edge(points_list, grid);
                } else {
                    paths = compute_paths(points_list, grid);
                }

                if (paths == null) {
                    continue;
                }

                if (GENERATE_INTERMEDIATE_RESULTS) {
                    draw(grid, paths, 0, false, iteration_location, width, scale);
                    //draw_paths(grid, paths, 0, false, iteration_location, width, scale);
                    //draw_flow(grid, paths, 0, false, iteration_location, width, scale);
                }

                ArrayList<double[][]> distances_for_paths = null;
                if (DISTANCE_METRIC.equals("BFS")) {
                    distances_for_paths = compute_bfs(grid, paths);
                } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                    distances_for_paths = compute_Dijkstra(grid, paths);
                } else if (DISTANCE_METRIC.equals("ANGULAR")) {
                    distances_for_paths = compute_angular_distance_precise(grid, paths);
                } else if (DISTANCE_METRIC.equals("ARC")) {
                    distances_for_paths = compute_arc_length(grid, paths);
                } else if (DISTANCE_METRIC.equals("ANGULAR_INTERSECTION")) {
                    distances_for_paths = compute_angular_distance_with_intersection(grid, paths);
                } else if (DISTANCE_METRIC.equals("ANGULAR_WITH_ARC_LENGTH")) {
                    distances_for_paths = compute_anguar_with_arc_length(grid, paths);
                } else if (DISTANCE_METRIC.equals("POLAR_SYSTEM")) {
                    distances_for_paths = compute_vertical_distances(grid, paths);
                }

                if (DRAW_DISTANCE_IMAGES) {
                    draw_distances(grid, paths, distances_for_paths, false, width, scale, 0, iteration_location);
                }

                // copies the initial state of grid into memory
                copy_height(grid, base_function);

                Tuple<Cell[][], ArrayList<Path>> result = iterate(grid, points_list, paths, distances_for_paths, NR_OF_ITERATIONS,
                        BASE_HEIGHT_TYPE, BASE_SCALE, true, false, width, scale, iteration_location, base_function);

                if (result == null) {
                    continue;
                }

                grid = result.first;

                paths = result.second;

                // draw final iteration
                draw(grid, paths, NR_OF_ITERATIONS, false, iteration_location, width, scale);
                draw_paths(grid, paths, NR_OF_ITERATIONS, false, iteration_location, width, scale);
                //draw_flow(grid, paths, NR_OF_ITERATIONS, false, iteration_location, width, scale);


                if (DRAW_DISTANCE_IMAGES) {
                    draw_distances(grid, paths, distances_for_paths, false, width, scale, NR_OF_ITERATIONS, iteration_location);
                }

                if (GENERATE_INTERMEDIATE_RESULTS) {
                    generate_gif(iteration_location, "global_height");
                }

                if (GENERATE_INTERMEDIATE_HEIGHT) {
                    generate_gif(iteration_location, "update_local_height");
                    generate_gif(iteration_location, "update_global_height");
                }

                log_file_writer.close();


            }
        }
    }

    public static void copy_height(Cell[][] grid, double[][] memory_grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                memory_grid[i][j] = grid[i][j].height;

            }
        }

    }

    public static void initialize_grid_height_to_edge(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {
                grid[j][i].height = NR_OF_COLUMNS - j;
            }
        }
    }

    public static void initialize_grid_height_to_edge_sqrt(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {
                grid[j][i].height = Math.sqrt(NR_OF_COLUMNS - j);
            }
        }
    }

    public static void initialize_grid_height_to_edge_squared(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {
                grid[j][i].height = Math.pow(NR_OF_COLUMNS - j, 2);
            }
        }
    }

    public static Cell[][] initialize_obstacles() {

        return null;
    }

    public static void compute_min_and_max_heights(Cell[][] grid) {

        for (int i = 0; i < NR_OF_ROWS; i++) {
            for (int j = 0; j < NR_OF_COLUMNS; j++) {

                if (grid[j][i].height < MIN_HEIGHT) {
                    MIN_HEIGHT = grid[j][i].height;
                }
                if (grid[j][i].height > MAX_HEIGHT) {
                    MAX_HEIGHT = grid[j][i].height;
                }

            }
        }
    }

    public static void initialize_parameters() {

        NR_OF_ROWS = 500;
        NR_OF_COLUMNS = 500;

        TARGET_NAME = "A";//"FL";
        INPUT_FILE_NAME = "./input/to_edge_4.csv";//"./input/1_s_20_t.csv";//"./input/1_s_8_t.csv";//"./input/USPos.csv";
        GIF_DELAY = 500; // 1000 - 1 FRAME PER SEC

        BASE_SCALE = 0.05;

        RESET_HEIGHTS = false;
        REMOVE_DIAGONAL_BIAS = false;

        DRAW_TEXT_DESCRIPTION = false;
        DRAW_PATHS = true;

        COLOR_MODE = "GRAY_SCALE";
        COLOR_MODE = "MULTIPLE_COLORS";
        //COLOR_MODE = "RED_BLUE";

        DRAW_DISTANCE_IMAGES = false;

        ARC_RADIUS = 20;

        BASE_HEIGHT_TYPE = "EUCLID";
        //BASE_HEIGHT_TYPE = "chebyshev";
        //BASE_HEIGHT_TYPE = "EUCLID_SQRT";
        //BASE_HEIGHT_TYPE = "EUCLID_SQUARED"; // previously known as default
        BASE_HEIGHT_TYPE = "TO_EDGE";
        //BASE_HEIGHT_TYPE = "TO_EDGE_SQUARED";
        //BASE_HEIGHT_TYPE = "TO_EDGE_SQRT";

        DISTANCE_METRIC = "DIJKSTRA";
        //DISTANCE_METRIC = "BFS";
        //DISTANCE_METRIC = "ANGULAR"; //  OLD!!!
        //DISTANCE_METRIC = "ARC";
        //DISTANCE_METRIC = "ANGULAR_INTERSECTION";
        //DISTANCE_METRIC = "ANGULAR_WITH_ARC_LENGTH";
        //DISTANCE_METRIC = "POLAR_SYSTEM";

        NR_OF_ITERATIONS = 10;

        WIDTHS = new double[]{100};
        SCALES = new double[]{10000};

        GENERATE_INTERMEDIATE_RESULTS = true;
        GENERATE_INTERMEDIATE_HEIGHT = true;

        EXPERIMENTAL_MODE = true;

        PATH_SCALING = true;
        SCALING_MODE = "WIDTHS";
        //SCALING_MODE = "FACTORS";

        MEMORY_MODE = false;
        MEMORY_DECAY_RATE = 0.66;

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
        fileWriter.write("ARC_RADIUS = " + ARC_RADIUS + "\n");
        fileWriter.write("RESET_HEIGHTS = " + RESET_HEIGHTS + "\n");
        fileWriter.write("REMOVE_DIAGONAL_BIAS = " + REMOVE_DIAGONAL_BIAS + "\n");
        fileWriter.write("TARGET_NAME = " + TARGET_NAME + "\n");
        fileWriter.write("PATH_SCALING = " + PATH_SCALING + "\n");
        fileWriter.write("MEMORY_MODE = " + MEMORY_MODE + "\n");
        fileWriter.write("MEMORY_DECAY_RATE = " + MEMORY_DECAY_RATE + "\n");

        fileWriter.close();
    }

    public static void generate_gif(String iteration_location, String type) throws IOException {

//        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/image_" + image_index + ".png"));
        File dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location + "\\"));

        String output_name = "";
        if (type.equals("global_height")) {
            output_name = "global_height";
        } else if (type.equals("update_local_height")) {
            output_name = "update_local_height";
        } else if (type.equals("update_global_height")) {
            output_name = "update_global_height";
        }

        String finalOutput_name1 = output_name;
        File[] files = dir.listFiles((dir1, name) -> (name.endsWith(".png") && name.startsWith(finalOutput_name1)));

        BufferedImage first = ImageIO.read(new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location + "\\" + output_name + "_0" + ".png")));
        ImageOutputStream output = new FileImageOutputStream(new File(currentWorkingPath.concat("\\" + storage_location_name + "\\" + iteration_location +
                "\\gif_" + output_name + ".gif")));

        GifSequenceWriter writer = new GifSequenceWriter(output, first.getType(), GIF_DELAY, true);
        writer.writeToSequence(first);


        String finalOutput_name = output_name;
        Arrays.sort(files, new Comparator<File>() {
            @Override
            public int compare(File f1, File f2) {
                String s1 = f1.getName().substring(finalOutput_name.length() + 1, f1.getName().indexOf("."));
                String s2 = f2.getName().substring(finalOutput_name.length() + 1, f2.getName().indexOf("."));
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

    public static void apply_height_update(Cell[][] grid, double[][] height_update) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[i][j].height = grid[i][j].height + height_update[i][j];

            }
        }

    }


    public static double[][] add_heights(double[][] previous_height, double[][] new_height) {

        double[][] result = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                // MEMORY_DECAY_RATE = 0.5; // ##hardcoded
                result[i][j] = previous_height[i][j] * MEMORY_DECAY_RATE + new_height[i][j];

            }
        }

        return result;
    }

    // adds the propagated height to the base and stores it in global grid
    public static Cell[][] add_to_base_function(Cell[][] grid, double[][] base, double[][] previous_height) {

        double[][] result = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[i][j].height = base[i][j] + previous_height[i][j];

            }
        }

        return grid;
    }

    public static Tuple<Cell[][], ArrayList<Path>> iterate(
            Cell[][] grid, ArrayList points_list, ArrayList paths,
            ArrayList distances_for_paths, int NR_OF_ITERATIONS,
            String BASE_HEIGHT_TYPE, double base_function_scale, boolean save_outputs,
            boolean show_intermediate_results, double width, double scale, String iteration_location,
            double[][] base_function)
            throws IOException {

        double[][] sum_of_previous_heights = new double[NR_OF_ROWS][NR_OF_COLUMNS];

        if (MEMORY_MODE) {
            // this is the first height update w.r.t. the base

            // take the old height (that is, the global grid height), compute height update U_1
            double[][] height_update = compute_height_update(grid, distances_for_paths, paths, 0, iteration_location); // this is U_1

            // here we apply the height update on the base function
            // the update is written inside global grid height
            apply_height_update(grid, height_update); // grid = B + U_1

            // sum of previous heights is now U_1
            sum_of_previous_heights = height_update;

        } else {
            update_height_overwrite(grid, distances_for_paths, width, scale, paths, 0, iteration_location);
        }

        compute_min_and_max_heights(grid);

        //adjust_height_min_distances(grid, distances_for_paths, width, scale, paths);

        // Iterate a number of times:
        for (int i = 1; i < NR_OF_ITERATIONS; i++) {
            System.out.println("iteration : " + i);
            log_file_writer.write("iteration : " + i + "\n");

            compute_flow(grid, i);

            if (EXPERIMENTAL_MODE) {
                paths = compute_paths_to_frame_edge(points_list, grid);
            } else {
                paths = compute_paths(points_list, grid);
            }

            if (paths == null) {
                return null;
            }

            if (save_outputs == true) {
                if (GENERATE_INTERMEDIATE_RESULTS) {
                    draw(grid, paths, i, show_intermediate_results, iteration_location, width, scale);
                    //draw_paths(grid, paths, i, show_intermediate_results, iteration_location, width, scale);
                    //draw_flow(grid, paths, i, show_intermediate_results, iteration_location, width, scale);

                }

                if (DRAW_DISTANCE_IMAGES) {
                    draw_distances(grid, paths, distances_for_paths, show_intermediate_results, width, scale, i, iteration_location);
                }
            }

            long startTime = System.currentTimeMillis();

            if (DISTANCE_METRIC.equals("BFS")) {
                distances_for_paths = compute_bfs(grid, paths);
            } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                distances_for_paths = compute_Dijkstra(grid, paths);
            } else if (DISTANCE_METRIC.equals("ANGULAR")) {
                distances_for_paths = compute_angular_distance_precise(grid, paths);
            } else if (DISTANCE_METRIC.equals("ARC")) {
                distances_for_paths = compute_arc_length(grid, paths);
            } else if (DISTANCE_METRIC.equals("ANGULAR_INTERSECTION")) {
                distances_for_paths = compute_angular_distance_with_intersection(grid, paths);
            } else if (DISTANCE_METRIC.equals("ANGULAR_WITH_ARC_LENGTH")) {
                distances_for_paths = compute_anguar_with_arc_length(grid, paths);
            } else if (DISTANCE_METRIC.equals("POLAR_SYSTEM")) {
                distances_for_paths = compute_vertical_distances(grid, paths);
            }

            long endTime = System.currentTimeMillis();

            Iterator path_it = distances_for_paths.iterator();

            ArrayList transposed_dist = new ArrayList();

            while (path_it.hasNext()) {

                DistanceForPathMatrix dist = (DistanceForPathMatrix) path_it.next();
                double[][] dist_matrix = dist.distance_matrix;

                transposeMatrix(dist_matrix);
                DistanceForPathMatrix dist_transposed = new DistanceForPathMatrix();
                dist_transposed.distance_matrix = dist_matrix;
                dist_transposed.path_id = dist.path_id;
                transposed_dist.add(dist_transposed);

            }

            distances_for_paths = transposed_dist;
            System.out.println("That took " + (endTime - startTime) + " milliseconds");
            log_file_writer.write("That took " + (endTime - startTime) + " milliseconds" + "\n");

            if (RESET_HEIGHTS == true) {
                if (BASE_HEIGHT_TYPE.equals("EUCLID")) {
                    initialize_grid_height_Euclidean_dist(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQUARED")) {
                    initialize_grid_height_Euclidean_squared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("chebyshev")) {
                    initialize_grid_height_chebyshev_distance(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQRT")) {
                    initialize_grid_height_Euclidean_dist_sqrt(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE")) {
                    initialize_grid_height_to_edge(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQUARED")) {
                    initialize_grid_height_to_edge_squared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQRT")) {
                    initialize_grid_height_to_edge_sqrt(grid);
                }
            }

            if (MEMORY_MODE) {
                // take the global grid height and compute the height update // here we took B + U_1 and computed U_2
                double[][] height_update = compute_height_update(grid, distances_for_paths, paths, i, iteration_location); // = U_2

                sum_of_previous_heights = add_heights(sum_of_previous_heights, height_update);

                // apply the height update on the global grid // here we should take B and add 1/2 U_1 + U_2
                //apply_height_update(grid, height_update); // this is not right

                grid = add_to_base_function(grid, base_function, sum_of_previous_heights);

            } else {
                update_height_overwrite(grid, distances_for_paths, width, scale, paths, i, iteration_location);
            }

            compute_min_and_max_heights(grid);

            //adjust_height_min_distances(grid, distances_for_paths, width, scale, paths);

        }

        compute_flow(grid, NR_OF_ITERATIONS);

        if (EXPERIMENTAL_MODE) {
            paths = compute_paths_to_frame_edge(points_list, grid);
        } else {
            paths = compute_paths(points_list, grid);
        }

        if (paths == null) {
            return null;
        }

        Tuple<Cell[][], ArrayList<Path>> tuple = new Tuple<Cell[][], ArrayList<Path>>(grid, paths);

        return tuple;

    }

    public static ArrayList compute_y_coordinates_for_columns(Cell[][] grid, ArrayList paths) {

        // TODO: potentially make this faster with binary search or something

        ArrayList y_coordinates_for_all_columns = new ArrayList();

        for (int col = 0; col < NR_OF_COLUMNS; col++) {

            ArrayList y_coordinates_for_column = new ArrayList();

            for (int row = 0; row < NR_OF_ROWS; row++) {

                Iterator path_iterator = paths.iterator();
                int intersection_counter = 0;

                while (path_iterator.hasNext()) {

                    Path path = (Path) path_iterator.next();

                    int path_id = path.id;

                    Iterator path_cell_iterator = path.cells.iterator();

                    while (path_cell_iterator.hasNext()) {

                        Cell path_cell = (Cell) path_cell_iterator.next();

                        if (grid[col][row] == path_cell) {
                            intersection_counter++;

                            // path cell is on the vertical cell
                            Pair<Integer, Integer> id_and_y_coord = new Pair<>(path_id, row);

                            y_coordinates_for_column.add(id_and_y_coord);
                            //path_indexes.add(path_index);
                        }
                    }
                }
            }

            //Map<Integer, List<Integer>> y_coord_hash_map = new HashMap<>();
            Map<Integer, List<Integer>> multiMap = new HashMap<Integer, List<Integer>>();

            for (int i = 0; i < y_coordinates_for_column.size(); i++) {

                Pair pair = (Pair) y_coordinates_for_column.get(i);

                List<Integer> list;
                if (multiMap.containsKey(pair.getKey())) {
                    list = multiMap.get(pair.getKey());
                    list.add((Integer) pair.getValue());
                } else {
                    list = new ArrayList<Integer>();
                    list.add((Integer) pair.getValue());
                    multiMap.put((Integer) pair.getKey(), list);
                }
            }

            for (int i = 0; i < multiMap.size(); i++) {
                if (multiMap.get(i).size() > 1) {

                    int sum = 0;
                    for (int j = 0; j < multiMap.get(i).size(); j++) {
                        sum = sum + multiMap.get(i).get(j);
                    }
                    int avg = sum / multiMap.get(i).size();

                    if (multiMap.get(i).size() > 3) {
                        //System.out.println();
                    }

                    Iterator y_iterator = y_coordinates_for_column.iterator();

                    ArrayList found = new ArrayList();
                    while (y_iterator.hasNext()) {

                        Pair pair = (Pair) y_iterator.next();

                        if ((Integer) pair.getKey() == i) {
                            found.add(pair);
                            //y_coordinates_for_column.remove(pair);
                            //e--;
                        }

                    }
                    y_coordinates_for_column.removeAll(found);

//                    int size = y_coordinates_for_column.size();
//                    for (int e = 0; e < size; e ++) {
//
//                        Pair pair = (Pair) y_coordinates_for_column.get(e);
//
//                        if((Integer) pair.getKey() == i) {
//                            y_coordinates_for_column.remove(e);
//                            //e--;
//                        }
//
//                    }

                    y_coordinates_for_column.add(i, new Pair<>(i, avg));

                }
            }

            y_coordinates_for_all_columns.add(y_coordinates_for_column);

        }

        return y_coordinates_for_all_columns;
    }

    public static double compute_scaled_path_factors(Cell[][] grid, ArrayList paths, ArrayList distances_for_cell, Cell cell) {

        ArrayList y_coordinates = new ArrayList();
        ArrayList path_indexes = new ArrayList();

        for (int i = 0; i < NR_OF_ROWS; i++) {

            Iterator path_iterator = paths.iterator();

            int path_index = 0;

            while (path_iterator.hasNext()) {

                ArrayList path = (ArrayList) path_iterator.next();

                Iterator path_cell_iterator = path.iterator();

                while (path_cell_iterator.hasNext()) {

                    Cell path_cell = (Cell) path_cell_iterator.next();

                    if (grid[cell.cell_col][i] == path_cell) {

                        // path cell is on the vertical cell
                        y_coordinates.add(i);
                        path_indexes.add(path_index);
                    }
                }
            }
        }

        ArrayList factors = new ArrayList();

        int index = -1;

        for (int i = 0; i < y_coordinates.size(); i++) {

            if (i == y_coordinates.size() - 1) {
                index = i;
                break;
            }

            if (cell.cell_row <= (int) y_coordinates.get(i + 1)) {
                index = i;
                break;
            }

            if (cell.cell_row <= (int) y_coordinates.get(i + 1) && cell.cell_row >= (int) y_coordinates.get(i)) {
                index = i;
                break;
            }
        }

        for (int i = 0; i < y_coordinates.size(); i++) {

            int distance_from_index = Math.abs(i - index);
            factors.add(1 / (Math.pow(2, distance_from_index)));

        }

        double sum = 0;
        for (int k = 0; k < distances_for_cell.size(); k++) {
            sum = sum + (double) factors.get(k) * gaussian((double) distances_for_cell.get(k), 0, HEIGHT_FUNCTION_WIDTH);

        }

        return sum;
    }

    public static double gaussian(double x, double mu, double sigma) {

        return (double) (1.0 / (Math.sqrt(2.0 * Math.PI) * sigma) * Math.exp(-Math.pow((x - mu) / sigma, 2.0) / 2));

    }

    public static void adjust_height_min_distances(Cell[][] grid, ArrayList distances_for_paths, double width, double scale, ArrayList paths) throws IOException {
        System.out.println("adjusting height");
        log_file_writer.write("adjusting height" + "\n");

        double[][] dist = (double[][]) distances_for_paths.get(0);

        double[][] min_distances = new double[NR_OF_ROWS][NR_OF_COLUMNS];

        for (int i = 0; i < NR_OF_ROWS; i++) {
            for (int j = 0; j < NR_OF_COLUMNS; j++) {

                Iterator dist_iter = distances_for_paths.iterator();

                double min_distance = (double) dist[i][j];

                while (dist_iter.hasNext()) {

                    double[][] distances_for_path = (double[][]) dist_iter.next();

                    if (distances_for_path[i][j] < min_distance) {

                        min_distance = distances_for_path[i][j];

                    }

                }
                min_distances[i][j] = min_distance;
            }
        }

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                Cell cell = grid[j][i];

                double new_height = gaussian((double) min_distances[cell.cell_row][cell.cell_col], 0, HEIGHT_FUNCTION_WIDTH);

                double height = ((-HEIGHT_FUNCTION_SCALE * new_height));
                grid[j][i].height = grid[j][i].height + ((height));

            }
        }

    }

    public static double[][] compute_height_update(Cell[][] grid, ArrayList distances_for_paths, ArrayList paths,
                                                   int iteration, String iteration_location) throws IOException {

        System.out.println("adjusting height");
        log_file_writer.write("adjusting height" + "\n");

        double[][] computed_height = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        ArrayList y_coordinates_for_columns = null;
        if (PATH_SCALING) {
            y_coordinates_for_columns = compute_y_coordinates_for_columns(grid, paths);
        }

        for (int row = 0; row < NR_OF_ROWS; row++) { // i is the row id
            //System.out.println();
            for (int col = 0; col < NR_OF_COLUMNS; col++) { //j is the column id

                if (PATH_SCALING) {

                    //SCALING_MODE = "FACTORS";

                    if (SCALING_MODE.equals("FACTORS")) {

                        factors_for_paths(grid, col, row, paths,
                                y_coordinates_for_columns, computed_height,
                                distances_for_paths);

                    } else if (SCALING_MODE.equals("WIDTHS")) {

                        widths_for_paths(grid, col, row, paths,
                                y_coordinates_for_columns, computed_height,
                                distances_for_paths);

                    }


                    // =================== PATH SCALING = FALSE
                    // ===================
                    // ===================
                    // ===================

                } else {


                    Cell cell = grid[col][row];

                    Iterator path_iterator = distances_for_paths.iterator();

                    ArrayList distances_for_cell = new ArrayList();

                    while (path_iterator.hasNext()) {

                        DistanceForPathMatrix distances = (DistanceForPathMatrix) path_iterator.next();
                        double[][] distances_matrix = distances.distance_matrix;

                        distances_for_cell.add(distances_matrix[cell.cell_row][cell.cell_col]);

                    }

                    double sum = 0;
                    for (int k = 0; k < distances_for_cell.size(); k++) {
                        sum = sum + gaussian((double) distances_for_cell.get(k), 0, HEIGHT_FUNCTION_WIDTH);
                    }
                    //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

                    long startTime = System.currentTimeMillis();

                    //double sum_2 = compute_scaled_path_factors(grid, paths, distances_for_cell, cell);

                    long endTime = System.currentTimeMillis();
                    //System.out.println("That took " + (endTime - startTime) + " milliseconds");

                    if (RESET_HEIGHTS == false) {

                        // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        //grid[col][row].height = grid[col][row].height + ((height));
                        computed_height[col][row] = height;

                    } else {
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        computed_height[col][row] = height;

                        //grid[col][row].height = grid[col][row].height + ((height));
                    }
                }
            }
        }


        if (GENERATE_INTERMEDIATE_RESULTS) {
            if (GENERATE_INTERMEDIATE_HEIGHT) {
                compute_min_and_max_heights(grid);
                draw_matrix(computed_height, paths, iteration, iteration_location, false);
                draw_matrix(computed_height, paths, iteration, iteration_location, true);

            }
        }

        return computed_height;
    }

    public static void widths_for_paths(Cell[][] grid, int col, int row, ArrayList paths,
                                        ArrayList y_coordinates_for_columns, double[][] computed_height,
                                        ArrayList distances_for_paths) {

        Cell cell = grid[col][row];

        // return the column of y_coordinate intersections of paths for the cell
        ArrayList y_coordinates_of_paths_for_cell = (ArrayList) y_coordinates_for_columns.get(cell.cell_col);

        HashMap abs_distances_for_cell = new HashMap();

        for (int k = 0; k < y_coordinates_of_paths_for_cell.size(); k++) {

            Pair pair = (Pair) y_coordinates_of_paths_for_cell.get(k);

            //Pair abs_d_pair = new Pair(pair.getKey(), Math.abs(row - (Integer) pair.getValue()));

            abs_distances_for_cell.put(pair.getKey(), Math.abs(row - (Integer) pair.getValue()));

        }

        //System.out.println();

        Iterator path_iterator = distances_for_paths.iterator();
        ArrayList distances_for_cell = new ArrayList();

        while (path_iterator.hasNext()) {

            DistanceForPathMatrix distanceForPathMatrix = (DistanceForPathMatrix) path_iterator.next();

            Pair distances_for_cell_and_id = new Pair(distanceForPathMatrix.path_id, distanceForPathMatrix.distance_matrix[cell.cell_row][cell.cell_col]);

            distances_for_cell.add(distances_for_cell_and_id);

        }

        double sum = 0;

        for (int k = 0; k < distances_for_cell.size(); k++) {

            Pair distances_for_cell_and_id = (Pair) distances_for_cell.get(k);

            double distance_for_path = (double) distances_for_cell_and_id.getValue();

            sum = sum + gaussian((double) distance_for_path, 0, 100 * ((int) abs_distances_for_cell.get(distances_for_cell_and_id.getKey()) + 1));

        }

        if (RESET_HEIGHTS == false) {

            // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));   // grid[col][row]
            grid[col][row].height = grid[col][row].height + ((height)); // updates horizontally
            computed_height[col][row] = height;

        } else {
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            computed_height[col][row] = height;

            grid[col][row].height = grid[col][row].height + ((height));
        }

    }

    public static void factors_for_paths(Cell[][] grid, int col, int row, ArrayList paths,
                                         ArrayList y_coordinates_for_columns, double[][] computed_height,
                                         ArrayList distances_for_paths) {

        Cell cell = grid[col][row];


        if (row == 5 && col == 5) {
            System.out.println();
        }

        if (row == 5 && col == 0) {
            System.out.println();
        }

        if (row == 9 && col == 0) {
            System.out.println();
        }

        if (row == 5 && col == 8) {
            System.out.println();
        }

        // grid[col][row]
        // TODO: potentially has to be sorted s.t. the path ids are preserved
        ArrayList y_coordinates_of_paths_for_cell = (ArrayList) y_coordinates_for_columns.get(cell.cell_col);

        //double[] factors = new double[paths.size()];

        Pair<Integer, Double>[] factors_for_ids = new Pair[paths.size()];

        //int index = -1;

        // TODO: this could be binary search (if sorted)
        // Check where the cell is w.r.t. the intersections of paths (y coordinates)
        for (int p = 0; p < y_coordinates_of_paths_for_cell.size(); p++) {

            Pair<Integer, Integer> pair = (Pair<Integer, Integer>) y_coordinates_of_paths_for_cell.get(p);

            if (p + 1 < y_coordinates_of_paths_for_cell.size()) {
                Pair<Integer, Integer> next_pair = (Pair<Integer, Integer>) y_coordinates_of_paths_for_cell.get(p + 1);

                if (row > (int) pair.getValue() && row < (int) next_pair.getValue()) {
                    // index is between two y coordinates
                    // the two y coordinates should have 1

                    factors_for_ids[p] = new Pair<>(pair.getKey(), 1.0);
                    factors_for_ids[p + 1] = new Pair<>(next_pair.getKey(), 1.0);
//
//                                factors[p] = 1;
//                                factors[p + 1] = 1;
                    break;
                }
            }

            if (p == (int) pair.getValue()) {
                // the index is on a path cell. We find all the paths that have the same y and set factors to 1
                //index = p;
//
//                            factors.add(p, 1);
//
                int counter = 0;
                Pair<Integer, Integer> next_pair = (Pair<Integer, Integer>) y_coordinates_of_paths_for_cell.get(p + counter);
                while (p + counter == next_pair.getValue()) {

                    factors_for_ids[p + counter] = new Pair<>(next_pair.getKey(), 1.0);
                    //factors[p + counter] = 1;

                    counter++;

                    if (p + counter == y_coordinates_of_paths_for_cell.size()) {
                        break;
                    }
                    next_pair = (Pair<Integer, Integer>) y_coordinates_of_paths_for_cell.get(p + counter);
                }
                break;

            }

            if (p < pair.getValue()) {
                //index = p;
                factors_for_ids[p] = new Pair<>(pair.getKey(), 1.0);
                //factors[p] = 1;
                break;
            }

        }

        int n = factors_for_ids.length;
        int first = -1;
        int last = -1;

        for (int k = 0; k < n; k++) {
//                        if () {
//
//                        }

            if (factors_for_ids[k] == null) {
                continue;
            }
            if (1.0 != factors_for_ids[k].getValue()) {
                continue;
            }

            if (first == -1)
                first = k;
            last = k;
        }

        for (int p = 0; p < y_coordinates_of_paths_for_cell.size(); p++) {

            int distance_from_index_first = Math.abs(p - first);
            int distance_from_index_last = Math.abs(p - last);

            double factor_first = 1 / (Math.pow(2, distance_from_index_first));
            double factor_last = 1 / (Math.pow(2, distance_from_index_last));

            //Pair new_pair = new Pair();
            Pair<Integer, Integer> pair = (Pair<Integer, Integer>) y_coordinates_of_paths_for_cell.get(p);

            //try {
            factors_for_ids[p] = new Pair<>(pair.getKey(), Math.max(factor_first, factor_last));
//                        } catch (Exception e) {
//                            System.out.println();
//                        }

        }
        // System.out.println();

        Iterator path_iterator = distances_for_paths.iterator();
        ArrayList distances_for_cell = new ArrayList();

        while (path_iterator.hasNext()) {

            DistanceForPathMatrix distanceForPathMatrix = (DistanceForPathMatrix) path_iterator.next();

            Pair distances_for_cell_and_id = new Pair(distanceForPathMatrix.path_id, distanceForPathMatrix.distance_matrix[cell.cell_row][cell.cell_col]);

            distances_for_cell.add(distances_for_cell_and_id);

        }

        double sum = 0;

        for (int k = 0; k < distances_for_cell.size(); k++) {

            //System.out.println();
            Pair distances_for_cell_and_id = (Pair) distances_for_cell.get(k);

            int path_id = (int) distances_for_cell_and_id.getKey();
            double distance_for_path = (double) distances_for_cell_and_id.getValue();
            double factor = 1;

            for (int b = 0; b < factors_for_ids.length; b++) {

                if (factors_for_ids[b] == null) {
                    factor = 1.0;
                    break;
                }

                if (factors_for_ids[b].getKey() == k) {
                    factor = factors_for_ids[b].getValue();
                    break;
                }

            }

            //double factor_for_path_id = factors_for_ids[];

            //factor = 1.0;
            //System.out.print("row : " + row + " col : " + col + " factor : " + factor + " | ");
            //factor = 1.0;
            sum = sum + factor * gaussian((double) distance_for_path, 0, HEIGHT_FUNCTION_WIDTH);

//                        if (factors[k] == 0) {
//                            factors[k] = 1;
//                        }
//
//                        sum = sum + factors[k] * gaussian((double) distances_for_cell.get(k), 0, HEIGHT_FUNCTION_WIDTH);
        }
        //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

        if (RESET_HEIGHTS == false) {

            // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));   // grid[col][row]
            grid[col][row].height = grid[col][row].height + ((height)); // updates horizontally
            computed_height[col][row] = height;

        } else {
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            computed_height[col][row] = height;

            grid[col][row].height = grid[col][row].height + ((height));
        }

    }

    public static void update_height_overwrite(Cell[][] grid, ArrayList distances_for_paths, double width, double scale, ArrayList paths, int iteration, String iteration_location) throws IOException {

        System.out.println("adjusting height");
        log_file_writer.write("adjusting height" + "\n");

        boolean verbose = false;

        double[][] computed_height = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        ArrayList y_coordinates_for_columns = null;
        if (PATH_SCALING) {
            y_coordinates_for_columns = compute_y_coordinates_for_columns(grid, paths);
        }

        for (int row = 0; row < NR_OF_ROWS; row++) { // i is the row id
            //System.out.println();
            for (int col = 0; col < NR_OF_COLUMNS; col++) { //j is the column id

                if (PATH_SCALING) {

                    if (SCALING_MODE.equals("FACTORS")) {

                        factors_for_paths(grid, col, row, paths,
                                y_coordinates_for_columns, computed_height,
                                distances_for_paths);

                    } else if (SCALING_MODE.equals("WIDTHS")) {

                        widths_for_paths(grid, col, row, paths,
                                y_coordinates_for_columns, computed_height,
                                distances_for_paths);

                    }


                    // =================== PATH SCALING = FALSE
                    // ===================
                    // ===================
                    // ===================

                } else {
                    Cell cell = grid[col][row];

                    Iterator path_iterator = distances_for_paths.iterator();

                    ArrayList distances_for_cell = new ArrayList();

                    while (path_iterator.hasNext()) {

                        DistanceForPathMatrix distances = (DistanceForPathMatrix) path_iterator.next();
                        double[][] distances_matrix = distances.distance_matrix;

                        distances_for_cell.add(distances_matrix[cell.cell_row][cell.cell_col]);

                    }

                    double sum = 0;
                    for (int k = 0; k < distances_for_cell.size(); k++) {
                        sum = sum + gaussian((double) distances_for_cell.get(k), 0, HEIGHT_FUNCTION_WIDTH);
                    }
                    //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

                    long startTime = System.currentTimeMillis();

                    //double sum_2 = compute_scaled_path_factors(grid, paths, distances_for_cell, cell);

                    long endTime = System.currentTimeMillis();
                    //System.out.println("That took " + (endTime - startTime) + " milliseconds");

                    if (RESET_HEIGHTS == false) {

                        // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        grid[col][row].height = grid[col][row].height + ((height));
                        computed_height[col][row] = height;

                    } else {
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        computed_height[col][row] = height;

                        grid[col][row].height = grid[col][row].height + ((height));
                    }
                }
            }
        }

        if (GENERATE_INTERMEDIATE_RESULTS) {
            if (GENERATE_INTERMEDIATE_HEIGHT) {
                compute_min_and_max_heights(grid);
                draw_matrix(computed_height, paths, iteration, iteration_location, false);
                draw_matrix(computed_height, paths, iteration, iteration_location, true);

            }
        }
    }

    public static void initialize_grid_height_Euclidean_squared(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (BASE_SCALE * (Math.pow(grid[j][i].cell_col - source_x, 2) + Math.pow(grid[j][i].cell_row - source_y, 2)));

            }
        }
    }

    public static void initialize_grid_height_Euclidean_dist_sqrt(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (BASE_SCALE * Math.sqrt(Math.sqrt(Math.pow(grid[j][i].cell_col - source_x, 2) + Math.pow(grid[j][i].cell_row - source_y, 2))));

            }
        }

    }

    public static void initialize_grid_height_Euclidean_dist_squared(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (BASE_SCALE * Math.pow(2, Math.sqrt(Math.pow(grid[j][i].cell_col - source_x, 2) + Math.pow(grid[j][i].cell_row - source_y, 2))));

            }
        }

    }

    public static void initialize_grid_height_chebyshev_distance(Cell[][] grid) {

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

        visited[source_cell.cell_row][source_cell.cell_col] = true;

        while (!queue.isEmpty()) {

            Cell cell = queue.peek();

            int x = cell.cell_col;
            int y = cell.cell_row;

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


    public static void initialize_grid_height_Euclidean_dist(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (BASE_SCALE * Math.sqrt(Math.pow(grid[j][i].cell_col - source_x, 2) + Math.pow(grid[j][i].cell_row - source_y, 2)));

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

    public static ArrayList compute_arc_length(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing arc length distance");
        log_file_writer.write("computing arc length distance" + "\n");

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

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - cell.cell_row, 2)));

                    if (radius > ARC_RADIUS) {

                        // Binary search finds the first cell for which the distance to this cell is >= than radius

//                        if (i == 389 && j == 0) {
//                            System.out.println();
//                        }

                        int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersection_cell = null;
                        //try {
                        intersection_cell = (Cell) path.get(index_of_cell);

//                        } catch (Exception e) {
//                            System.out.println();
//                        }

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - intersection_cell.cell_row, 2)));

                        // TODO: check index out of bounds
                        if (index_of_cell + 1 < path.size() && index_of_cell - 1 > 0) {

                            Cell next_cell = (Cell) path.get(index_of_cell + 1);
                            Cell previous_cell = (Cell) path.get(index_of_cell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(source_cell.cell_col - next_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - next_cell.cell_row, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(source_cell.cell_col - previous_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - previous_cell.cell_row, 2)));

                            Tuple<Double, Double> intersection_point = null;
                            // if radius is between intersection_cell and intersection_cell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        previous_cell.cell_col, previous_cell.cell_row);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersection_cell and intersection_cell + 1 we consider these two cells

                                // here compute the intersection

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        next_cell.cell_col, next_cell.cell_row);

                            } else if (radius == dist) {
                                intersection_point = new Tuple<Double, Double>((double) intersection_cell.cell_col, (double) intersection_cell.cell_row);
                            }

                            if (intersection_point == null) {
                                //System.out.println();
                            }

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(cell.cell_row - intersection_point.second, 2)));

                            double distance_from_source_to_intersection = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(source_cell.cell_row - intersection_point.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distance_from_source_to_intersection, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * distance_from_source_to_intersection));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length;

                        } else {

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_cell.cell_col, 2) +
                                    Math.pow(cell.cell_row - intersection_cell.cell_row, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * radius));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length;

                        }

                    } else {

                        int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersection_cell = (Cell) path.get(index_of_cell);

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - intersection_cell.cell_row, 2)));

                        // TODO: check index out of bounds
                        if (index_of_cell + 1 < path.size() && index_of_cell - 1 > 0) {

                            Cell next_cell = (Cell) path.get(index_of_cell + 1);
                            Cell previous_cell = (Cell) path.get(index_of_cell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(source_cell.cell_col - next_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - next_cell.cell_row, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(source_cell.cell_col - previous_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - previous_cell.cell_row, 2)));

                            Tuple<Double, Double> intersection_point = null;
                            // if radius is between intersection_cell and intersection_cell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        previous_cell.cell_col, previous_cell.cell_row);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersection_cell and intersection_cell + 1 we consider these two cells

                                // here compute the intersection

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        next_cell.cell_col, next_cell.cell_row);

                            } else if (radius == dist) {
                                intersection_point = new Tuple<Double, Double>((double) intersection_cell.cell_col, (double) intersection_cell.cell_row);
                            }

                            if (intersection_point == null) {
                                //System.out.println();
                            }

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(cell.cell_row - intersection_point.second, 2)));

                            double distance_from_source_to_intersection = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(source_cell.cell_row - intersection_point.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distance_from_source_to_intersection, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * distance_from_source_to_intersection));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length;

                        } else {

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_cell.cell_col, 2) +
                                    Math.pow(cell.cell_row - intersection_cell.cell_row, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * radius));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length;
                        }
                    }
                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_col][cell.cell_row] = 0.0;

            }

            double[][] transposed_matrix = transposeMatrix(distances);

            distances_for_paths.add(transposed_matrix);

        }
        return distances_for_paths;
    }

    public static ArrayList compute_angular_with_Dijkstra(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular distance with arc length ");
        log_file_writer.write("computing angular distance with arc length " + "\n");

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

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - cell.cell_row, 2)));

                    if (radius > ARC_RADIUS) {

                        // Binary search finds the first cell for which the distance to this cell is >= than radius
                        int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersection_cell = (Cell) path.get(index_of_cell);

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - intersection_cell.cell_row, 2)));

                        // TODO: check index out of bounds
                        if (index_of_cell + 1 < path.size() && index_of_cell - 1 > 0) {

                            Cell next_cell = (Cell) path.get(index_of_cell + 1);
                            Cell previous_cell = (Cell) path.get(index_of_cell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(source_cell.cell_col - next_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - next_cell.cell_row, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(source_cell.cell_col - previous_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - previous_cell.cell_row, 2)));

                            Tuple<Double, Double> intersection_point = null;
                            // if radius is between intersection_cell and intersection_cell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        previous_cell.cell_col, previous_cell.cell_row);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersection_cell and intersection_cell + 1 we consider these two cells

                                // here compute the intersection

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        next_cell.cell_col, next_cell.cell_row);

                            } else if (radius == dist) {
                                intersection_point = new Tuple<Double, Double>((double) intersection_cell.cell_col, (double) intersection_cell.cell_row);
                            }

                            if (intersection_point == null) {
                                //System.out.println();
                            }

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(cell.cell_row - intersection_point.second, 2)));

                            double distance_from_source_to_intersection = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(source_cell.cell_row - intersection_point.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distance_from_source_to_intersection, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * distance_from_source_to_intersection));

                            distances[i][j] = angle;

                        } else {

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_cell.cell_col, 2) +
                                    Math.pow(cell.cell_row - intersection_cell.cell_row, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * radius));

                            distances[i][j] = angle;

                        }

                    } else {

                        int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersection_cell = (Cell) path.get(index_of_cell);

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - intersection_cell.cell_row, 2)));

                        // TODO: check index out of bounds
                        if (index_of_cell + 1 < path.size() && index_of_cell - 1 > 0) {

                            Cell next_cell = (Cell) path.get(index_of_cell + 1);
                            Cell previous_cell = (Cell) path.get(index_of_cell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(source_cell.cell_col - next_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - next_cell.cell_row, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(source_cell.cell_col - previous_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - previous_cell.cell_row, 2)));

                            Tuple<Double, Double> intersection_point = null;
                            // if radius is between intersection_cell and intersection_cell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        previous_cell.cell_col, previous_cell.cell_row);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersection_cell and intersection_cell + 1 we consider these two cells

                                // here compute the intersection

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        next_cell.cell_col, next_cell.cell_row);

                            } else if (radius == dist) {
                                intersection_point = new Tuple<Double, Double>((double) intersection_cell.cell_col, (double) intersection_cell.cell_row);
                            }

                            if (intersection_point == null) {
                                // System.out.println();
                            }

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(cell.cell_row - intersection_point.second, 2)));

                            double distance_from_source_to_intersection = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(source_cell.cell_row - intersection_point.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distance_from_source_to_intersection, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * distance_from_source_to_intersection));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length / radius;

                        } else {

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_cell.cell_col, 2) +
                                    Math.pow(cell.cell_row - intersection_cell.cell_row, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * radius));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length / radius;
                        }

                    }
                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_col][cell.cell_row] = 0.0;

            }

            distances_for_paths.add(transposeMatrix(distances));
        }
        return distances_for_paths;

    }

    public static ArrayList compute_anguar_with_arc_length(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular distance with arc length ");
        log_file_writer.write("computing angular distance with arc length " + "\n");

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

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - cell.cell_row, 2)));

                    if (radius > ARC_RADIUS) {

                        // Binary search finds the first cell for which the distance to this cell is >= than radius
                        int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersection_cell = (Cell) path.get(index_of_cell);

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - intersection_cell.cell_row, 2)));

                        // TODO: check index out of bounds
                        if (index_of_cell + 1 < path.size() && index_of_cell - 1 > 0) {

                            Cell next_cell = (Cell) path.get(index_of_cell + 1);
                            Cell previous_cell = (Cell) path.get(index_of_cell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(source_cell.cell_col - next_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - next_cell.cell_row, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(source_cell.cell_col - previous_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - previous_cell.cell_row, 2)));

                            Tuple<Double, Double> intersection_point = null;
                            // if radius is between intersection_cell and intersection_cell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        previous_cell.cell_col, previous_cell.cell_row);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersection_cell and intersection_cell + 1 we consider these two cells

                                // here compute the intersection

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        next_cell.cell_col, next_cell.cell_row);

                            } else if (radius == dist) {
                                intersection_point = new Tuple<Double, Double>((double) intersection_cell.cell_col, (double) intersection_cell.cell_row);
                            }

                            if (intersection_point == null) {
                                // System.out.println();
                            }

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(cell.cell_row - intersection_point.second, 2)));

                            double distance_from_source_to_intersection = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(source_cell.cell_row - intersection_point.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distance_from_source_to_intersection, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * distance_from_source_to_intersection));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = angle;

                        } else {

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_cell.cell_col, 2) +
                                    Math.pow(cell.cell_row - intersection_cell.cell_row, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * radius));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = angle;

                        }

                    } else {

                        int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersection_cell = (Cell) path.get(index_of_cell);

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - intersection_cell.cell_row, 2)));

                        // TODO: check index out of bounds
                        if (index_of_cell + 1 < path.size() && index_of_cell - 1 > 0) {

                            Cell next_cell = (Cell) path.get(index_of_cell + 1);
                            Cell previous_cell = (Cell) path.get(index_of_cell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(source_cell.cell_col - next_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - next_cell.cell_row, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(source_cell.cell_col - previous_cell.cell_col, 2) +
                                    Math.pow(source_cell.cell_row - previous_cell.cell_row, 2)));

                            Tuple<Double, Double> intersection_point = null;
                            // if radius is between intersection_cell and intersection_cell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        previous_cell.cell_col, previous_cell.cell_row);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersection_cell and intersection_cell + 1 we consider these two cells

                                // here compute the intersection

                                intersection_point = compute_intersection_of_circle_and_line_segment(
                                        source_x, source_y, radius,
                                        intersection_cell.cell_col, intersection_cell.cell_row,
                                        next_cell.cell_col, next_cell.cell_row);

                            } else if (radius == dist) {
                                intersection_point = new Tuple<Double, Double>((double) intersection_cell.cell_col, (double) intersection_cell.cell_row);
                            }

                            if (intersection_point == null) {
                                //System.out.println();
                            }

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(cell.cell_row - intersection_point.second, 2)));

                            double distance_from_source_to_intersection = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_point.first, 2) +
                                    Math.pow(source_cell.cell_row - intersection_point.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distance_from_source_to_intersection, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * distance_from_source_to_intersection));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length / ARC_RADIUS;

                        } else {

                            double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_cell.cell_col, 2) +
                                    Math.pow(cell.cell_row - intersection_cell.cell_row, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                    (2.0 * radius * radius));

                            double arc_length = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arc_length / ARC_RADIUS;
                        }
                    }
                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_col][cell.cell_row] = 0.0;

            }

            double[][] transposed_matrix = transposeMatrix(distances);

            distances_for_paths.add(transposed_matrix);

//            for (int i = 0; i < NR_OF_COLUMNS; i++) {
//                for (int j = 0; j < NR_OF_ROWS; j++) {
//
//                    Cell cell = grid[i][j];
//
//                    double radius = (Math.sqrt(Math.pow(source_cell.cell_x - cell.cell_x, 2) +
//                            Math.pow(source_cell.cell_y - cell.cell_y, 2)));
//
//                    if (radius < ARC_RADIUS) {
//                        //System.out.print("r : " + r + " c : " + c + " " + ANSI_YELLOW + first[r][c] + ANSI_RESET + " ");
//
//                        System.out.print(" " + ANSI_YELLOW + distances[i][j] + ANSI_RESET);
//                    } else {
//                        System.out.print(" " + distances[i][j]);
//                    }
//
//                }
//                System.out.println();
//            }

        }
        return distances_for_paths;

    }

    public static ArrayList compute_angular_distance_with_intersection(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular intersections distance");
        log_file_writer.write("computing angular intersections distance" + "\n");

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

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - cell.cell_row, 2)));

                    int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                    Cell intersection_cell = (Cell) path.get(index_of_cell);

                    double dist = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - intersection_cell.cell_row, 2)));

                    // TODO: check index out of bounds
                    if (index_of_cell + 1 < path.size() && index_of_cell - 1 > 0) {
                        Cell next_cell = null;
                        //try {
                        next_cell = (Cell) path.get(index_of_cell + 1);

                        //} catch (Exception e) {
                        //    System.out.println();
                        // }
                        Cell previous_cell = (Cell) path.get(index_of_cell - 1);

                        double dist_1 = (Math.sqrt(Math.pow(source_cell.cell_col - next_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - next_cell.cell_row, 2)));

                        double dist_2 = (Math.sqrt(Math.pow(source_cell.cell_col - previous_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - previous_cell.cell_row, 2)));

                        Tuple<Double, Double> intersection_point = null;
                        // if radius is between intersection_cell and intersection_cell - 1, we consider these two cells
                        if (radius > dist_2 && radius < dist) {

                            // here compute the intersection point

                            intersection_point = compute_intersection_of_circle_and_line_segment(
                                    source_x, source_y, radius,
                                    intersection_cell.cell_col, intersection_cell.cell_row,
                                    previous_cell.cell_col, previous_cell.cell_row);


                        } else if (radius > dist && radius < dist_1) {
                            // if radius is between intersection_cell and intersection_cell + 1 we consider these two cells

                            // here compute the intersection

                            intersection_point = compute_intersection_of_circle_and_line_segment(
                                    source_x, source_y, radius,
                                    intersection_cell.cell_col, intersection_cell.cell_row,
                                    next_cell.cell_col, next_cell.cell_row);

                        } else if (radius == dist) {
                            intersection_point = new Tuple<Double, Double>((double) intersection_cell.cell_col, (double) intersection_cell.cell_row);
                        }

                        if (intersection_point == null) {
                            // System.out.println();
                        }

                        double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_point.first, 2) +
                                Math.pow(cell.cell_row - intersection_point.second, 2)));

                        double distance_from_source_to_intersection = (Math.sqrt(Math.pow(source_cell.cell_col - intersection_point.first, 2) +
                                Math.pow(source_cell.cell_row - intersection_point.second, 2)));

                        double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distance_from_source_to_intersection, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                (2.0 * radius * distance_from_source_to_intersection));

                        distances[i][j] = angle;

                    } else {

                        double distance_from_cell_to_intersection = (Math.sqrt(Math.pow(cell.cell_col - intersection_cell.cell_col, 2) +
                                Math.pow(cell.cell_row - intersection_cell.cell_row, 2)));

                        double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distance_from_cell_to_intersection, 2)) /
                                (2.0 * radius * radius));

                        distances[i][j] = angle;

                    }

                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_col][cell.cell_row] = 0.0;

            }

            //distances_for_paths.add(distances);
            distances_for_paths.add(transposeMatrix(distances));
        }
        return distances_for_paths;


    }


    public static Tuple<Double, Double> compute_intersection_of_circle_and_line_segment(double circlex, double circley, double radius,
                                                                                        double x1, double y1,
                                                                                        double x2, double y2) {

        //Calculate change in x and y for the segment
        double deltax = x2 - x1;
        double deltay = y2 - y1;

        //Set up our quadratic formula
        double a = deltax * deltax + deltay * deltay;
        double b = 2 * (deltax * (x1 - circlex) + deltay * (y1 - circley));
        double c = (x1 - circlex) * (x1 - circlex) + (y1 - circley) * (y1 - circley) - radius * radius;

        //Check if there is a negative in the discriminant
        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0)
            return null;

        //Try both +- in the quadratic formula
        double quad1 = (-b + Math.sqrt(discriminant)) / (2 * a);
        double quad2 = (-b - Math.sqrt(discriminant)) / (2 * a);

        //If the result is between 0 and 1, there is an intersection
        if (quad1 >= 0 && quad1 <= 1) {
            //System.out.println("quad1 : " + quad1 + " quad2 : " + quad2);

            double x = x1 + quad1 * (x2 - x1);
            double y = y1 + quad1 * (y2 - y1);

            //System.out.println("x = " + x + " y = " + y);

            return new Tuple<>(x, y);
        } else if (quad2 >= 0 && quad2 <= 1) {
            //System.out.println("quad1 : " + quad1 + " quad2 : " + quad2);

            double x = x1 + quad2 * (x2 - x1);
            double y = y1 + quad2 * (y2 - y1);

            //System.out.println("x = " + x + " y = " + y);

            return new Tuple<>(x, y);
        }
        return null;
    }

    public static ArrayList compute_vertical_distances(Cell[][] grid, ArrayList<Path> paths) throws IOException {
        System.out.println("computing vertical distances ");
        log_file_writer.write("computing vertical distances " + "\n");

        Iterator path_iterator = paths.iterator();

        ArrayList distances_for_paths = new ArrayList();

        int path_counter = 0;

        // for all paths
        while (path_iterator.hasNext()) {

            Path path = (Path) path_iterator.next();

            ArrayList path_cells = path.cells;

            Collections.reverse(path_cells);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    // inefficient method:

                    Iterator path_cell_iterator = path_cells.iterator();

                    while (path_cell_iterator.hasNext()) {

                        Cell path_cell = (Cell) path_cell_iterator.next();

                        if (path_cell.cell_col == cell.cell_col) {

                            double distance = Math.abs(path_cell.cell_row - cell.cell_row);
                            distances[i][j] = distance;
                        }

                    }

                }
            }

            Iterator path_cell_iter = path_cells.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_col][cell.cell_row] = 0.0;

            }

            double[][] transposed_distances_matrix = transposeMatrix(distances);
            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distance_matrix = transposed_distances_matrix;
            distanceForPathMatrix.path_id = path_counter;
            path_counter++;
            distances_for_paths.add(distanceForPathMatrix);
        }
        return distances_for_paths;

    }

    public static ArrayList compute_Dijkstra(Cell[][] grid, ArrayList<Path> paths) throws IOException {

        System.out.println("Computing Dijkstra");
        log_file_writer.write("Computing Dijkstra" + "\n");

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

            Path path = (Path) path_iterator.next();

            // store results for a single path here
            //double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distance_matrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];
            distanceForPathMatrix.path_id = path.id;

            // Initialize distances with max distance values
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {
                    distanceForPathMatrix.distance_matrix[j][i] = Integer.MAX_VALUE;
                    grid[j][i].distance = Integer.MAX_VALUE;
                }
            }

            PriorityQueue<Cell> Q = new PriorityQueue<Cell>(NR_OF_COLUMNS * NR_OF_ROWS,
                    new distanceComparator());

            // Loop over each cell on a path:
            Iterator cell_iterator = path.cells.iterator();
            while (cell_iterator.hasNext()) {

                // A cell on a path:
                Cell cell = (Cell) cell_iterator.next();

                //cell.distance = 0.0;
                distanceForPathMatrix.distance_matrix[cell.cell_col][cell.cell_row] = 0.0;
                cell.distance = 0.0;

                Q.add(cell);

            }

            while (!Q.isEmpty()) {

                //Tuple<Cell, Double> nd = Q.poll();
                Cell cell = Q.poll();

                int x = cell.cell_col;
                int y = cell.cell_row;

                for (int i = 0; i < 16; i++) {

                    int adj_x = x + dCol[i];
                    int adj_y = y + dRow[i];
                    double weight = weights[i];

                    if (is_inside_grid(adj_x, adj_y)) {

                        if (distanceForPathMatrix.distance_matrix[adj_x][adj_y] > distanceForPathMatrix.distance_matrix[x][y] + weight) {

                            // If Cell is already been reached once,
                            // remove it from priority queue
                            if (distanceForPathMatrix.distance_matrix[adj_x][adj_y] != Integer.MAX_VALUE) {
                                Cell adj = grid[adj_x][adj_y];//new Cell(rows, cols, dist[rows][cols]);
                                adj.distance = distanceForPathMatrix.distance_matrix[adj_x][adj_y];
                                Q.remove(adj);

                            }

                            // Insert cell with updated distance
                            distanceForPathMatrix.distance_matrix[adj_x][adj_y] = ((distanceForPathMatrix.distance_matrix[cell.cell_col][cell.cell_row] + weight));

                            grid[adj_x][adj_y].distance = distanceForPathMatrix.distance_matrix[adj_x][adj_y];
                            Q.add(grid[adj_x][adj_y]); //new Cell(rows, cols, dist[rows][cols]));

                        }
                    }
                }
            }
            distanceForPathMatrix.distance_matrix = transposeMatrix(distanceForPathMatrix.distance_matrix);
            distances_for_paths.add(distanceForPathMatrix);
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

    public static ArrayList compute_angular_distance_precise(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular distance");
        log_file_writer.write("computing angular distance" + "\n");

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

                int cell_x = cell.cell_col;
                int next_cell_x = next_cell.cell_col;

                int cell_y = cell.cell_row;
                int next_cell_y = next_cell.cell_row;

//                double distance_x = Math.sqrt(Math.pow(next_cell.cell_x - cell.cell_x, 2));
//                double distance_y = Math.sqrt(Math.pow(next_cell.cell_y - cell.cell_y, 2));

                double abs_x = Math.abs(cell_x - next_cell_x);
                double abs_y = Math.abs(cell_y - next_cell_y);

                double splits = 1;

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
            last_cell.cell_x = last_path_cell.cell_col;
            last_cell.cell_y = last_path_cell.cell_row;

            finer_path.add(last_cell);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - cell.cell_row, 2)));

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
                    double dist = (Math.sqrt(Math.pow(intersection_cell.cell_x - source_cell.cell_col, 2) +
                            Math.pow(intersection_cell.cell_y - source_cell.cell_row, 2)));

                    double dist_2 = (Math.sqrt(Math.pow(intersection_cell.cell_x - cell.cell_col, 2) +
                            Math.pow(intersection_cell.cell_y - cell.cell_row, 2)));

                    double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(dist_2, 2)) /
                            (2.0 * radius * dist));

                    distances[i][j] = angle;
                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_col][cell.cell_row] = 0.0;

            }

            distances_for_paths.add(transposeMatrix(distances));
        }
        return distances_for_paths;

    }

    public static ArrayList compute_angular_distance(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular distance");
        log_file_writer.write("computing angular distance" + "\n");

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

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - cell.cell_row, 2)));

                    int index_of_cell = binary_search_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                    Cell intersection_cell = (Cell) path.get(index_of_cell);

                    // Here we can either use radius as dist or the actual distance. Matter of precision.
                    double dist = (Math.sqrt(Math.pow(intersection_cell.cell_col - source_cell.cell_col, 2) +
                            Math.pow(intersection_cell.cell_row - source_cell.cell_row, 2)));

                    double dist_2 = (Math.sqrt(Math.pow(intersection_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(intersection_cell.cell_row - cell.cell_row, 2)));

                    double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(dist_2, 2)) /
                            (2.0 * radius * dist));

                    distances[i][j] = angle;
                }
            }

            Iterator path_cell_iter = path.iterator();

            while (path_cell_iter.hasNext()) {

                Cell cell = (Cell) path_cell_iter.next();

                distances[cell.cell_col][cell.cell_row] = 0.0;

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

            double dist_to_mid = (Math.sqrt(Math.pow(source_cell.cell_col - mid_cell.cell_x, 2) +
                    Math.pow(source_cell.cell_row - mid_cell.cell_y, 2)));

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

            double dist_to_mid = (Math.sqrt(Math.pow(source_cell.cell_col - mid_cell.cell_col, 2) +
                    Math.pow(source_cell.cell_row - mid_cell.cell_row, 2)));

            if (dist_to_mid == distance_target) {
                return mid;
            }

            if (dist_to_mid <= distance_target) {
                start = mid + 1;
            } else {
                ans = mid;

                if (mid - 1 > 0) {
                    end = mid - 1;
                } else {
                    return mid;
                }
            }
        }
        return end;
    }

    public static int binarySearch(ArrayList<Cell> path, int l, int r, double dist) {

        if (r >= l) {
            int mid = l + (r - l) / 2;

            Cell mid_cell = path.get(mid);

            double dist_to_mid = (Math.sqrt(Math.pow(source_cell.cell_col - mid_cell.cell_col, 2) +
                    Math.pow(source_cell.cell_row - mid_cell.cell_row, 2)));

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

                    double radius = (Math.sqrt(Math.pow(source_cell.cell_col - cell.cell_col, 2) +
                            Math.pow(source_cell.cell_row - cell.cell_row, 2)));

                    Iterator cell_iterator = path.iterator();

                    while (cell_iterator.hasNext()) {

                        Cell path_cell = (Cell) cell_iterator.next();

                        double dist = (Math.sqrt(Math.pow(source_cell.cell_col - path_cell.cell_col, 2) +
                                Math.pow(source_cell.cell_row - path_cell.cell_row, 2)));
                        if (dist <= radius) {
                            // consider next cell on path
                            continue;
                        } else {
                            // this is the cell
                            double dist_2 = (Math.sqrt(Math.pow(cell.cell_col - path_cell.cell_col, 2) +
                                    Math.pow(cell.cell_row - path_cell.cell_row, 2)));

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

    public static ArrayList compute_bfs(Cell[][] grid, ArrayList<Path> paths) throws IOException {

        System.out.println("computing bfs");
        log_file_writer.write("computing bfs" + "\n");

        Iterator path_iterator = paths.iterator();

        // Direction vectors
        int dRow[] = {-1, 0, 1, 0};
        int dCol[] = {0, 1, 0, -1};

        // dw = weight of the edge

        ArrayList distances_for_paths = new ArrayList();

        while (path_iterator.hasNext()) {


            boolean[][] visited = new boolean[NR_OF_COLUMNS][NR_OF_ROWS];

            //double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distance_matrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    distanceForPathMatrix.distance_matrix[i][j] = 0.0f;

                }
            }

            Path path = (Path) path_iterator.next();

            distanceForPathMatrix.path_id = path.id;

            Queue<Cell> queue = new LinkedList<>();

            Iterator cell_iterator = path.cells.iterator();

            while (cell_iterator.hasNext()) {

                Cell cell = (Cell) cell_iterator.next();

                // add all cells of a path to queue
                queue.add(cell);

                visited[cell.cell_row][cell.cell_col] = true;

            }

            while (!queue.isEmpty()) {

                Cell cell = queue.peek();

                int x = cell.cell_col;
                int y = cell.cell_row;

                queue.remove();

                for (int i = 0; i < 4; i++) {

                    int adj_x = x + dCol[i];
                    int adj_y = y + dRow[i];

                    if (isValid(visited, adj_y, adj_x)) {

                        queue.add(grid[adj_x][adj_y]);
                        visited[adj_y][adj_x] = true;

                        distanceForPathMatrix.distance_matrix[adj_y][adj_x] = distanceForPathMatrix.distance_matrix[y][x] + 1;
                    }
                }
            }

            distances_for_paths.add(distanceForPathMatrix);
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

                        double distance = (double) (Math.sqrt(Math.pow(path_cell.cell_col - current_cell.cell_col, 2) + Math.pow(path_cell.cell_row - current_cell.cell_row, 2)));

                        distances.add(distance);

                    }

                    int min_distance_index = distances.indexOf(Collections.min(distances));

                    Cell min_distance_cell = (Cell) path.get(min_distance_index);

                }
            }
        }
    }

    public static void draw_distances(Cell[][] grid, ArrayList paths, ArrayList distances_for_paths,
                                      boolean show_intermediate_results, double width, double scale, int image_index,
                                      String iteration_location) throws
            IOException {

//        jframe = new JFrame("panel");
//        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);
//
//        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
//                BufferedImage.TYPE_INT_ARGB);
//
//        double[][] dist = (double[][]) distances_for_paths.get(1);
//
//        double max_dist = dist[0][0];
//        double min_dist = 0;
//
//        for (int i = 0; i < NR_OF_COLUMNS; i++) {
//
//            for (int j = 0; j < NR_OF_ROWS; j++) {
//
//                if (dist[i][j] > max_dist) {
//                    max_dist = dist[i][j];
//                }
//            }
//        }
//
//        for (int i = 0; i < NR_OF_COLUMNS; i++) {
//            for (int j = 0; j < NR_OF_ROWS; j++) {
//
//                int c = (int) (dist[i][j] * 255.0 / max_dist);
//
//
//                if (c < 0) {
//                    image.setRGB(i, j, new Color(0, 0, -c).getRGB());
//                } else {
//                    image.setRGB(i, j, new Color(c, 0, 0).getRGB());
//                }
//
//
////                if (!grid[i][j].title.equals("")) {
////
////                    image.setRGB(i, j, new Color(0, 255, 0).getRGB());
////                }
//            }
//        }
//
//        Iterator iter = paths.iterator();
//
//        while (iter.hasNext()) {
//
//            ArrayList path = (ArrayList) iter.next();
//
//            Iterator cell_iter = path.iterator();
//
//            while (cell_iter.hasNext()) {
//
//                Cell cell = (Cell) cell_iter.next();
//
//                image.setRGB((int) cell.cell_x, (int) cell.cell_y, new Color(255, 255, 255).getRGB());
//
//            }
//        }
//
//        JPanel pane = new JPanel() {
//            @Override
//            protected void paintComponent(Graphics g) {
//                super.paintComponent(g);
//                g.drawImage(image, 0, 0, null);
//            }
//        };
//
//        jframe.add(pane);
//
//        jframe.setVisible(true);
//        jframe.show();
//
//        ImageIO.write(image, "png", new File("image.png"));


        //////////////////// ++++++++++++++++++++


        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        double[][] dist = (double[][]) distances_for_paths.get(0);

        double max_dist = dist[0][0];
        double min_dist = 0;

        double[][] min_distances = new double[NR_OF_ROWS][NR_OF_COLUMNS];

        for (int i = 0; i < NR_OF_ROWS; i++) {
            for (int j = 0; j < NR_OF_COLUMNS; j++) {

                Iterator dist_iter = distances_for_paths.iterator();

                double min_distance = (double) dist[i][j];

                while (dist_iter.hasNext()) {

                    double[][] distances_for_path = (double[][]) dist_iter.next();

                    if (distances_for_path[i][j] < min_distance) {

                        min_distance = distances_for_path[i][j];

                    }

                }
                min_distances[i][j] = min_distance;

            }
        }

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                if (min_distances[i][j] > max_dist) {
                    max_dist = min_distances[i][j];
                }
            }
        }

        System.out.println("min dist: " + min_dist + " max dist : " + max_dist);
        log_file_writer.write("min dist: " + min_dist + " max dist : " + max_dist + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                //int c = (int) (dist[i][j] * 255.0 / max_dist);

                float value = (float) ((min_distances[j][i] - min_dist) / (max_dist - min_dist));

                Color color = getValueBetweenTwoFixedColors(value);

                image.setRGB(i, j, color.getRGB());

                if (!grid[i][j].title.equals("")) {

                    image.setRGB(i, j, new Color(0, 255, 0).getRGB());
                }
            }
        }

        if (DRAW_PATHS) {

            Iterator iter = paths.iterator();

            while (iter.hasNext()) {

                Path path = (Path) iter.next();

                Iterator cell_iter = path.cells.iterator();

                while (cell_iter.hasNext()) {

                    Cell cell = (Cell) cell_iter.next();

                    image.setRGB((int) cell.cell_col, (int) cell.cell_row, new Color(255, 255, 255).getRGB());

                }
            }
        }

        // draw string in image:

        if (DRAW_TEXT_DESCRIPTION) {
            Font f = new Font(Font.MONOSPACED, Font.PLAIN, 20);
            String s = "width: " + width + " scale: " + scale + " i: " + image_index;
            Graphics g = image.getGraphics();
            g.setColor(Color.BLUE);
            g.setFont(f);
            FontMetrics fm = g.getFontMetrics();
            int x = image.getWidth() - fm.stringWidth(s) - 5;
            int y = fm.getHeight();
            g.drawString(s, x, y);
            g.dispose();
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

        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/image_distances_" + image_index + ".png"));
        file.mkdirs();
        ImageIO.write(image, "png", file);

    }

    public static Color getValueBetweenTwoFixedColors(float value) {

        int aR;
        int aG;
        int aB;
        int bR;
        int bG;
        int bB;

        if (COLOR_MODE == "GRAY_SCALE") {
            aR = 0;
            aG = 0;
            aB = 0;

            bR = 255;
            bG = 255;
            bB = 255;
        } else {
            aR = 0;
            aG = 0;
            aB = 255;  // RGB for our 1st color (blue in this case).

            bR = 255;
            bG = 0;
            bB = 0;    // RGB for our 2nd color (red in this case).
        }

        int red = (int) ((float) (bR - aR) * value + aR);      // Evaluated as -255*value + 255.
        int green = (int) ((float) (bG - aG) * value + aG);      // Evaluates as 0.
        int blue = (int) ((float) (bB - aB) * value + aB);      // Evaluates as 255*value + 0.

        Color color = null;
        //try {
        color = new Color(red, green, blue);
        // } catch (Exception e) {
        //     System.out.println();
        // }

        return color;
    }

    public static Color getHeatMapColor(float value) {

        // TODO: fix, low priority
        int NUM_COLORS = 5;

        float color[][] = {{0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}};
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

        Color result_color = new Color(red, green, blue);

        return result_color;
    }

    public static void draw_matrix(double[][] matrix, ArrayList paths, int image_index, String iteration_location, boolean relative_to_total) throws IOException {
        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        double max_height = matrix[0][0];
        double min_hieght = matrix[0][0];

        if (relative_to_total) {
            max_height = MAX_HEIGHT;
            min_hieght = MIN_HEIGHT;

            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    if (matrix[i][j] < min_hieght) {
                        min_hieght = matrix[i][j];
                    }

                    if (matrix[i][j] > max_height) {
                        max_height = matrix[i][j];
                    }

                }
            }

        } else {
            for (int i = 0; i < NR_OF_COLUMNS; i++) {

                for (int j = 0; j < NR_OF_ROWS; j++) {

                    if (matrix[i][j] > max_height) {
                        max_height = matrix[i][j];
                    }
                    if (matrix[i][j] < min_hieght) {
                        min_hieght = matrix[i][j];
                    }
                }
            }
        }

        System.out.println("min height update : " + min_hieght + " max height update : " + max_height);
        log_file_writer.write("min height update : " + min_hieght + " max height update : " + max_height + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                float value = (float) ((matrix[i][j] - min_hieght) / (max_height - min_hieght));

                Color color = null;
                if (COLOR_MODE == "GRAY_SCALE") {
                    color = getValueBetweenTwoFixedColors(value);

                } else if (COLOR_MODE == "MULTIPLE_COLORS") {
                    float minHue = 210f / 255;//210f / 255; //corresponds to green
                    float maxHue = 0; //corresponds to red
                    float hue = value * maxHue + (1 - value) * minHue;
                    color = new Color(Color.HSBtoRGB(hue, 1f, 1f)); //getHeatMapColor(value);
                } else if (COLOR_MODE == "RED_BLUE") {
                    color = getValueBetweenTwoFixedColors(value);
                }

                // try {
                image.setRGB(i, j, color.getRGB());

                //  } catch (Exception e) {
                //      System.out.println();
                //  }

            }
        }

//        if (DRAW_PATHS) {
//
//            Iterator iter = paths.iterator();
//
//            while (iter.hasNext()) {
//
//                Path path = (Path) iter.next();
//
//                Iterator cell_iter = path.iterator();
//
//                while (cell_iter.hasNext()) {
//
//                    Cell cell = (Cell) cell_iter.next();
//
//                    image.setRGB((int) cell.cell_x, (int) cell.cell_y, new Color(255, 255, 255).getRGB());
//
//                }
//            }
//        }
        File file;
        if (relative_to_total) {
            file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/update_global_height_" + image_index + ".png"));
        } else {
            file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/update_local_height_" + image_index + ".png"));
        }
        file.mkdirs();
        ImageIO.write(image, "png", file);
    }

    public static void draw_flow(Cell[][] grid, ArrayList paths, int image_index, boolean show_intermediate_results, String iteration_location, double width, double scale)
            throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                int flow = grid[i][j].flow_direction;
                Color color = null;

                if (flow == 1) {
                    color = Color.RED;
                } else if (flow == 2) {
                    color = Color.PINK;
                } else if (flow == 4) {
                    color = Color.ORANGE;
                } else if (flow == 8) {
                    color = Color.YELLOW;
                } else if (flow == 16) {
                    color = Color.GREEN;
                } else if (flow == 32) {
                    color = Color.CYAN;
                } else if (flow == 64) {
                    color = Color.BLUE;
                } else if (flow == 128) {
                    color = Color.MAGENTA;
                }

                image.setRGB(i, j, color.getRGB());

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

        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/paths_" + image_index + ".png"));
        file.mkdirs();
        ImageIO.write(image, "png", file);

    }

    public static void draw_paths(Cell[][] grid, ArrayList paths, int image_index, boolean show_intermediate_results, String iteration_location, double width, double scale)
            throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                Color color = Color.BLACK;

                image.setRGB(i, j, color.getRGB());

            }
        }

        if (DRAW_PATHS) {

            Iterator iter = paths.iterator();

            while (iter.hasNext()) {

                Path path = (Path) iter.next();

                Iterator cell_iter = path.cells.iterator();

                while (cell_iter.hasNext()) {

                    Cell cell = (Cell) cell_iter.next();

                    image.setRGB((int) cell.cell_col, (int) cell.cell_row, new Color(255, 255, 255).getRGB());

                }
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

        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/paths_" + image_index + ".png"));
        file.mkdirs();
        ImageIO.write(image, "png", file);

    }

    public static void draw(Cell[][] grid, ArrayList paths, int image_index, boolean show_intermediate_results, String iteration_location, double width, double scale)
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

        System.out.println("min height : " + min_hieght + " max height : " + max_height);
        log_file_writer.write("min height : " + min_hieght + " max height : " + max_height + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                float value = (float) ((grid[i][j].height - min_hieght) / (max_height - min_hieght));

                Color color = null;
                if (COLOR_MODE == "GRAY_SCALE") {
                    color = getValueBetweenTwoFixedColors(value);

                } else if (COLOR_MODE == "MULTIPLE_COLORS") {
                    float minHue = 210f / 255;//210f / 255; //corresponds to green
                    float maxHue = 0; //corresponds to red
                    float hue = value * maxHue + (1 - value) * minHue;
                    color = new Color(Color.HSBtoRGB(hue, 1f, 1f)); //getHeatMapColor(value);
                } else if (COLOR_MODE == "RED_BLUE") {
                    color = getValueBetweenTwoFixedColors(value);
                }

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


//                if (!grid[i][j].title.equals("")) {
//
//                    image.setRGB(i, j, new Color(0, 255, 0).getRGB());
//                }

            }
        }

        if (DRAW_PATHS) {

            Iterator iter = paths.iterator();

            while (iter.hasNext()) {

                Path path = (Path) iter.next();

                Iterator cell_iter = path.cells.iterator();

                while (cell_iter.hasNext()) {

                    Cell cell = (Cell) cell_iter.next();

                    image.setRGB((int) cell.cell_col, (int) cell.cell_row, new Color(255, 255, 255).getRGB());

                }
            }
        }


        // draw string in image:

        if (DRAW_TEXT_DESCRIPTION) {
            Font f = new Font(Font.MONOSPACED, Font.PLAIN, 20);
            String s = "width: " + width + " scale: " + scale + " i: " + image_index;
            Graphics g = image.getGraphics();
            g.setColor(Color.BLUE);
            g.setFont(f);
            FontMetrics fm = g.getFontMetrics();
            int x = image.getWidth() - fm.stringWidth(s) - 5;
            int y = fm.getHeight();
            g.drawString(s, x, y);
            g.dispose();
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

        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/" + iteration_location + "/global_height_" + image_index + ".png"));
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

    public static ArrayList compute_paths_to_frame_edge(ArrayList points_list, Cell[][] grid) throws IOException {
        System.out.println("computing paths from fram to edge");
        log_file_writer.write("computing paths from fram to edge" + "\n");

        Iterator it = points_list.iterator();

        // ArrayList<ArrayList<Cell>> paths = new ArrayList();
        ArrayList<Path> paths = new ArrayList();

        int path_id_counter = 0;
        while (it.hasNext()) {

            Point point = (Point) it.next();

            if (point.name.equals(TARGET_NAME)) {
                continue;
            }

            Cell grid_cell = grid[point.grid_x][point.grid_y];

            Cell current_cell = grid_cell;

            ArrayList<Cell> path = new ArrayList();

            path.add(current_cell);

            int counter = 0;

            while (!(current_cell.title.equals("right_edge"))) {

                if (counter > 4 * (NR_OF_COLUMNS + NR_OF_ROWS)) {
                    System.out.println("something went wrong");
                    log_file_writer.write("something went wrong" + "\n");
                    return null;
                }

                int node_x = (int) current_cell.cell_col;
                int node_y = (int) current_cell.cell_row;

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

                counter++;

            }

            Path new_path = new Path();
            new_path.cells = path;
            new_path.id = path_id_counter;
            paths.add(new_path);
            path_id_counter++;
        }
        return paths;
    }

    public static ArrayList<Path> compute_paths(ArrayList points_list, Cell[][] grid) throws IOException {

        System.out.println("computing paths");
        log_file_writer.write("computing paths" + "\n");

        Iterator it = points_list.iterator();

        ArrayList<Path> paths = new ArrayList();

        int id_counter = 0;
        while (it.hasNext()) {

            Point point = (Point) it.next();

            if (point.name.equals(TARGET_NAME)) {
                continue;
            }

            Cell grid_cell = grid[point.grid_x][point.grid_y];

            Cell current_cell = grid_cell;

            Path path = new Path();
            path.cells = new ArrayList();
            path.id = id_counter;
            //ArrayList<Cell> path = new ArrayList();

            path.cells.add(current_cell);

            int counter = 0;

            while (!(current_cell.title.equals(TARGET_NAME))) {

                if (counter > 4 * (NR_OF_COLUMNS + NR_OF_ROWS)) {
                    System.out.println("something went wrong");
                    log_file_writer.write("something went wrong" + "\n");

                    return null;
                }

                int node_x = (int) current_cell.cell_col;
                int node_y = (int) current_cell.cell_row;

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

                path.cells.add(current_cell);

                counter++;

            }
            id_counter++;

            paths.add(path);

        }
        return paths;
    }

    public static void compute_flow(Cell[][] grid, int iteration) throws IOException {
        System.out.println("computing flow");
        log_file_writer.write("computing flow" + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                int x = (int) grid[i][j].cell_col;
                int y = (int) grid[i][j].cell_row;

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

                    double distance = 0.0;

                    if (neighbor == left || neighbor == right || neighbor == top || neighbor == bottom) {

                        distance = 1.0;

                    } else if (neighbor == top_left || neighbor == top_right || neighbor == bottom_left || neighbor == bottom_right) {

                        double dist_from_source_to_cell = Math.sqrt(Math.pow(grid[j][i].cell_col - source_x, 2) + Math.pow(grid[j][i].cell_row - source_y, 2));

                        if (true) {//dist_from_source_to_cell < 100) {

                            double difference = Math.sqrt(2) - 1.0;
                            double step = difference / NR_OF_ITERATIONS;

                            if (REMOVE_DIAGONAL_BIAS) {

                                //System.out.println("bias removed");
                                distance = Math.sqrt(2) - step * iteration;
                                ;//Math.sqrt(2) - step * iteration;//1.0;
                            } else {
                                //System.out.println("bias not removed");
                                distance = (double) Math.sqrt(2);
                            }

                        } else {
                            //System.out.println("bias not removed ");

                            distance = (double) Math.sqrt(2);
                        }
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
                grid[point.grid_x][point.grid_y].height = -10000;

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

                cell.cell_col = i;
                cell.cell_row = j;
                cell.flow_direction = 0;
                cell.height = 0;

                cell.title = "";
                grid[i][j] = cell;
                cell.is_obstacle = false;

                // i = columns
                // j = rows
                if (i == nr_of_rows - 1 || j == nr_of_columns - 1 || i == 0 || j == 0) {
                    cell.title = "edge";
                }

                if (i == nr_of_rows - 1) {
                    cell.title = "right_edge";
                }

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
