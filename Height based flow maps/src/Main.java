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

    public static int sourceX;
    public static int sourceY;

    public static String TARGET_NAME;

    static JFrame jframe;


    private static double zoomFactor = 50;
    private static double prevZoomFactor = 50;
    private static boolean zoomer;


    public static DateTimeFormatter dtf = DateTimeFormatter.ofPattern("MM_dd_HH_mm_ss");
    public static LocalDateTime now = LocalDateTime.now();
    public static String storageLocationName = "";
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

    public static Cell sourceCell;

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
    public static double previousSum = 0.0;
    public static double newSum = 0.0;
    public static boolean DRAW_DISTANCE_IMAGES = false;
    public static boolean HEIGHT_DECAY_ENABLED = false;
    public static double MIN_HEIGHT = Integer.MAX_VALUE;
    public static double MAX_HEIGHT = Integer.MIN_VALUE;
    public static boolean GENERATE_INTERMEDIATE_RESULTS = true;
    public static boolean GENERATE_INTERMEDIATE_HEIGHT = true;
    public static boolean HORIZONTAL_FLOW_MODE = false;
    public static FileWriter logFileWriter;
    public static String COLOR_MODE = "";
    public static boolean PATH_SCALING = false;
    public static boolean MEMORY_MODE = false;
    public static double MEMORY_DECAY_RATE;
    public static String SCALING_MODE;

    public static boolean FLOW_ACCUMULATION;

    public static boolean CIRCULAR_MODE;

    public static String OBSTACLES;

    public static int OBSTACLE_WIDTH;
    public static int OBSTACLE_SCALE;
    public static int OBSTACLE_RATE;

    public static double [][] HEIGHT_UPDATE;

    public static void main(String[] args) throws IOException {

        initializeParameters();

        String[] parts = INPUT_FILE_NAME.split("/");
        String testCaseName = parts[2].split("\\.")[0];

        storageLocationName = dtf.format(now);

        storageLocationName = storageLocationName.concat("_" + testCaseName);
        storageLocationName = storageLocationName.concat("_" + DISTANCE_METRIC + "_" + BASE_HEIGHT_TYPE);

        currentWorkingPath = System.getProperty("user.dir").concat("\\experiments\\");

        File dir = new File(currentWorkingPath.concat("\\" + storageLocationName + "\\"));
        dir.mkdir();

        for (int i = 0; i < WIDTHS.length; i++) {

            double width = WIDTHS[i];
            HEIGHT_FUNCTION_WIDTH = width;

            for (int j = 0; j < SCALES.length; j++) {
                double scale = SCALES[j];
                HEIGHT_FUNCTION_SCALE = scale;

                String iterationLocation = ("w_" + width + "_s_" + scale);

                dir = new File(currentWorkingPath.concat("\\" + storageLocationName + "\\" + iterationLocation + "\\"));

                dir.mkdir();

                logFileWriter = new FileWriter(currentWorkingPath.concat("\\" + storageLocationName + "\\" + iterationLocation) + "\\log.txt");

                writeOutputConfiguration(iterationLocation);

                Cell[][] grid = initializeGrid(NR_OF_ROWS, NR_OF_COLUMNS);
                double[][] baseFunction = new double[NR_OF_ROWS][NR_OF_COLUMNS];

                ArrayList<String> items = readInput();

                ArrayList<Point> pointsList = processInput(items);

                Bounds bounds = obtainBounds(pointsList);

                computeCellForPoint(bounds, pointsList);

                initializePointsInGrid(grid, pointsList);

                if (BASE_HEIGHT_TYPE.equals("EUCLID")) {
                    initializeGridHeight_EuclideanDist(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQUARED")) {
                    initializeGridHeight_EuclideanSquared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("chebyshev")) {
                    initializeGridHeightChebyshevDistance(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQRT")) {
                    initializeGridHeight_EuclideanDistSqrt(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE")) {
                    initializeGridHeightToEdge(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQUARED")) {
                    initializeGridHeightToEdgeSquared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQRT")) {
                    initializeGridHeightToEdgeSqrt(grid);
                }

                double[][] heightUpdateObstacles;
                computeMinAndMaxHeights(grid);
                if (OBSTACLES.equals("STATIC")) {
                    heightUpdateObstacles = computeObstaclesStatic(grid, pointsList);
                    drawObstacles(heightUpdateObstacles, new ArrayList(), 0, iterationLocation);

                } else if (OBSTACLES.equals("PROGRESSIVE")) {
                    heightUpdateObstacles = initializeObstaclesProgressive(grid, pointsList, 0);
                    drawObstacles(heightUpdateObstacles, new ArrayList(), 0, iterationLocation);
                }
                computeMinAndMaxHeights(grid);
                //computeMinAndMaxHeights(grid);

                computeFlow(grid, 0);
                //computeFlowAccumulation(grid);

                ArrayList<Path> paths;

                if (HORIZONTAL_FLOW_MODE) {
                    paths = computePathsToFrameEdge(pointsList, grid);
                } else {
                    paths = computePaths(pointsList, grid);
                }

                if (paths == null) {
                    continue;
                }

                if (GENERATE_INTERMEDIATE_RESULTS) {
                    draw(grid, paths, 0, false, iterationLocation, width, scale, pointsList);
                    //drawPaths(grid, paths, 0, false, iterationLocation, width, scale);
                    //drawFlow(grid, paths, 0, false, iterationLocation, width, scale);
                    //drawFlowAccumulation(grid, paths, 0, false, iterationLocation, width, scale);
                }

                ArrayList<DistanceForPathMatrix> distancesForPaths = null;
                if (DISTANCE_METRIC.equals("BFS")) {
                    distancesForPaths = computeBfs(grid, paths);
                } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                    distancesForPaths = compute_Dijkstra(grid, paths);
                } else if (DISTANCE_METRIC.equals("ARC")) {
                    distancesForPaths = computeArcLength(grid, paths);
                } else if (DISTANCE_METRIC.equals("ANGULAR_INTERSECTION")) {
                    distancesForPaths = computeAngularDistanceWithIntersection(grid, paths);
                } else if (DISTANCE_METRIC.equals("ANGULAR_WITH_ARC_LENGTH")) {
                    distancesForPaths = computeAnguarWithArcLength(grid, paths);
                } else if (DISTANCE_METRIC.equals("POLAR_SYSTEM")) {
                    distancesForPaths = computeVerticalDistances(grid, paths);
                } else if (DISTANCE_METRIC.equals("BRUTE_FORCE")) {
                    distancesForPaths = computeDistancesBruteForce(grid, paths);
                }

//                ArrayList<DistanceForPathMatrix> distancesForPaths1 = computeDistancesBruteForce(grid, paths);
//                ArrayList<DistanceForPathMatrix> distancesForPaths2 = compute_Dijkstra(grid, paths);
//
//                DistanceForPathMatrix distanceForPathMatrix1 = distancesForPaths1.get(0);
//                DistanceForPathMatrix distanceForPathMatrix2 = distancesForPaths2.get(0);
//
//                double sum1 = 0.0;
//
//                for (int c = 0 ; c < NR_OF_COLUMNS; c++) {
//
//                    for (int r = 0 ; r < NR_OF_ROWS ; r ++) {
//
//                        double difference = distanceForPathMatrix1.distanceMatrix[c][r] - distanceForPathMatrix2.distanceMatrix[c][r];
//
//                        sum1 = sum1 + Math.abs(difference);
//
//                        //System.out.print(difference + " ");
//                    }
//                   // System.out.println();
//                }
//
//                System.out.println("Sum of differences DIJKSTRA-BRUTE: " + sum1);
//
//                ArrayList<DistanceForPathMatrix> distancesForPaths3 = computeBfs(grid, paths);
//
//                DistanceForPathMatrix distanceForPathMatrix3 = distancesForPaths3.get(0);
//
//                double sum2 = 0.0;
//
//                for (int c = 0 ; c < NR_OF_COLUMNS; c++) {
//
//                    for (int r = 0 ; r < NR_OF_ROWS ; r ++) {
//
//                        double difference = distanceForPathMatrix1.distanceMatrix[c][r] - distanceForPathMatrix3.distanceMatrix[c][r];
//
//                        sum2 = sum2 + Math.abs(difference);
//
//                        //System.out.print(difference + " ");
//                    }
//                   // System.out.println();
//                }
//
//                System.out.println("Sum of differences BFS-BRUTE: " + sum2);
//
//                System.out.println("BFS is " + sum2 / sum1 + " times less precise than Dijkstra");

                if (DRAW_DISTANCE_IMAGES) {
                    drawDistances(grid, paths, distancesForPaths, false, width, scale, 0, iterationLocation);
                }

                // copies the initial state of grid into memory

                DistanceForPathMatrix d = (DistanceForPathMatrix) distancesForPaths.get(0);

//                for (int r = 0; r < NR_OF_ROWS; r++) {
//
//                    for (int c = 0; c < NR_OF_COLUMNS; c++) {
//
//                        double value = Math.round(d.distanceMatrix[r][c] * 10.0) / 10.0;
//                        System.out.print(value + " & ");
//
//                    }
//                    System.out.println(" \\\\ ");
//                }

                copyHeight(grid, baseFunction);

                Tuple<Cell[][], ArrayList<Path>> result = iterate(grid, pointsList, paths, distancesForPaths, NR_OF_ITERATIONS,
                        BASE_HEIGHT_TYPE, BASE_SCALE, true, false, width, scale, iterationLocation, baseFunction);

                if (result == null) {
                    continue;
                }

                grid = result.first;

                paths = result.second;

                // draw final iteration
                draw(grid, paths, NR_OF_ITERATIONS, false, iterationLocation, width, scale, pointsList);
                drawPaths(grid, paths, NR_OF_ITERATIONS, false, iterationLocation, width, scale);
                //drawFlow(grid, paths, NR_OF_ITERATIONS, false, iterationLocation, width, scale);
                //drawFlowAccumulation(grid, paths, NR_OF_ITERATIONS, false, iterationLocation, width, scale);


                if (DRAW_DISTANCE_IMAGES) {
                    drawDistances(grid, paths, distancesForPaths, false, width, scale, NR_OF_ITERATIONS, iterationLocation);
                }

                if (GENERATE_INTERMEDIATE_RESULTS) {
                    //generateGif(iterationLocation, "\"/height_w_\" + width + \"_s_\" + scale + \"_\"");
                }

                if (GENERATE_INTERMEDIATE_HEIGHT) {
                  //  generateGif(iterationLocation, "updateLocalHeight");
                 //   generateGif(iterationLocation, "updateGlobalHeight");
                }

                logFileWriter.close();


            }
        }
    }

    public static void copyHeight(Cell[][] grid, double[][] memoryGrid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                memoryGrid[i][j] = grid[i][j].height;

            }
        }

    }

    public static void initializeGridHeightToEdge(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {
                grid[j][i].height = NR_OF_COLUMNS - j;
            }
        }
    }

    public static void initializeGridHeightToEdgeSqrt(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {
                grid[j][i].height = Math.sqrt(NR_OF_COLUMNS - j);
            }
        }
    }

    public static void initializeGridHeightToEdgeSquared(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {
                grid[j][i].height = Math.pow(NR_OF_COLUMNS - j, 2);
            }
        }
    }

    public static double[][] computeObstaclesStatic(Cell[][] grid, ArrayList<Point> pointsList) {

        double[][] heightUpdate = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int n = 0; n < pointsList.size(); n++) {

            if (!(pointsList.get(n).name.equals(TARGET_NAME) || pointsList.get(n).name.equals("S"))) {

                Point point = pointsList.get(n);

                for (int col = 0; col < NR_OF_COLUMNS; col++) {
                    for (int row = 0; row < NR_OF_ROWS; row++) {

                        Cell cell = grid[col][row];

                        double distance = Math.sqrt(Math.pow(grid[col][row].cellCol - point.gridCol, 2) + Math.pow(grid[col][row].cellRow - point.gridRow, 2));

                        double width = OBSTACLE_WIDTH; // 20
                        double height = OBSTACLE_SCALE; //500

                        double obstacleHeight = gaussian(distance, 0, width);

                        double update = height * obstacleHeight;

                        heightUpdate[col][row] = heightUpdate[col][row] + update;

                        grid[col][row].height = grid[col][row].height + update;

                    }
                }
            }
        }
        return heightUpdate;
    }

    public static double[][] initializeObstaclesProgressive(Cell[][] grid, ArrayList<Point> pointsList, int iteration) {

        //double iterationFactor = (iteration) / NR_OF_ITERATIONS;

        //double fX = 1.0 - 1.0 / (1.0 + iteration);//(iteration + 1) / NR_OF_ITERATIONS;
        //double fMaxX = 1.0 - 1.0 / (1.0 + NR_OF_ITERATIONS);//(iteration + 1) / NR_OF_ITERATIONS;
        double gX = iteration * iteration * 0.01;//fX / fMaxX;

        System.out.println("Obstacle factor : " + gX);
        //iterationFactor = 1;

        double[][] heightUpdate = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int n = 0; n < pointsList.size(); n++) {

            if (!(pointsList.get(n).name.equals(TARGET_NAME) || pointsList.get(n).name.equals("S"))) {

                Point point = pointsList.get(n);

                for (int col = 0; col < NR_OF_COLUMNS; col++) {
                    for (int row = 0; row < NR_OF_ROWS; row++) {

                        Cell cell = grid[col][row];

                        double distance = Math.sqrt(Math.pow(grid[col][row].cellCol - point.gridCol, 2) + Math.pow(grid[col][row].cellRow - point.gridRow, 2));

                        double width = OBSTACLE_WIDTH;  // 20
                        double height = OBSTACLE_SCALE; // 500

                        double obstacleHeight = gaussian(distance, 0, width);

                        double update = gX * height * obstacleHeight;

                        heightUpdate[col][row] = heightUpdate[col][row] + update;

                        grid[col][row].height = grid[col][row].height + update;

                    }
                }
            }
        }
        return heightUpdate;
    }

    public static void computeMinAndMaxHeights(Cell[][] grid) {

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

    public static void initializeParameters() {

        NR_OF_ROWS = 500;
        NR_OF_COLUMNS = 500;
        HEIGHT_UPDATE = new double[NR_OF_COLUMNS][NR_OF_ROWS];
        TARGET_NAME = "A";//"FL";
        INPUT_FILE_NAME = "./input/obstacles.csv";//"./input/1S_20T.csv";//"./input/1S_8T.csv";//"./input/USPos.csv";
        GIF_DELAY = 500; // 1000 - 1 FRAME PER SEC

        BASE_SCALE = 0.05;//1.0;//0.05;//1.0;//0.05;

        RESET_HEIGHTS = true;
        REMOVE_DIAGONAL_BIAS = false;

        DRAW_TEXT_DESCRIPTION = false;
        DRAW_PATHS = true;

        //COLOR_MODE = "GRAY_SCALE";
        COLOR_MODE = "MULTIPLE_COLORS";
        //COLOR_MODE = "RED_BLUE";

        DRAW_DISTANCE_IMAGES = false;

        ARC_RADIUS = 20;

        BASE_HEIGHT_TYPE = "EUCLID";
        //BASE_HEIGHT_TYPE = "chebyshev";
        //BASE_HEIGHT_TYPE = "EUCLID_SQRT";
        //BASE_HEIGHT_TYPE = "EUCLID_SQUARED"; // previously known as default
        //BASE_HEIGHT_TYPE = "TO_EDGE";
        //BASE_HEIGHT_TYPE = "TO_EDGE_SQUARED";
        //BASE_HEIGHT_TYPE = "TO_EDGE_SQRT";

        DISTANCE_METRIC = "DIJKSTRA";
        //DISTANCE_METRIC = "BRUTE_FORCE";
        //DISTANCE_METRIC = "BFS";
        //DISTANCE_METRIC = "ANGULAR"; //  OLD!!!
        //DISTANCE_METRIC = "ARC";
        //DISTANCE_METRIC = "ANGULAR_INTERSECTION";
        //DISTANCE_METRIC = "ANGULAR_WITH_ARC_LENGTH";
        //DISTANCE_METRIC = "POLAR_SYSTEM";

        NR_OF_ITERATIONS = 100;

        WIDTHS = new double[]{20};
        SCALES = new double[]{10};

        GENERATE_INTERMEDIATE_RESULTS = true;
        GENERATE_INTERMEDIATE_HEIGHT = true;

        HORIZONTAL_FLOW_MODE = false;

        PATH_SCALING = false;
        SCALING_MODE = "WIDTHS";
        //SCALING_MODE = "OVERLAPS";

        FLOW_ACCUMULATION = false;

        MEMORY_MODE = false;
        MEMORY_DECAY_RATE = 0.66;

        CIRCULAR_MODE = false;

        OBSTACLES = "STATIC";
        //OBSTACLES = "PROGRESSIVE";
        //OBSTACLES = "";

        OBSTACLE_RATE = 1;
        OBSTACLE_WIDTH = 10;
        OBSTACLE_SCALE = 100; //10

    }

    public static void writeOutputConfiguration(String iterationLocation) throws IOException {
        ////        File file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/image_" + imageIndex + ".png"));
        FileWriter fileWriter = new FileWriter(currentWorkingPath.concat("\\" + storageLocationName + "\\" + iterationLocation) + "\\config.txt");
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
        fileWriter.write("OBSTACLE_WIDTH = " + OBSTACLE_WIDTH + "\n");
        fileWriter.write("OBSTACLE_SCALE = " + OBSTACLE_SCALE + "\n");
        fileWriter.write("OBSTACLE_RATE = " + OBSTACLE_RATE + "\n");
        fileWriter.write("OBSTACLES = " + OBSTACLES + "\n");
        fileWriter.close();
    }

    public static void generateGif(String iterationLocation, String type) throws IOException {

//        File file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/image_" + imageIndex + ".png"));
        File dir = new File(currentWorkingPath.concat("\\" + storageLocationName + "\\" + iterationLocation + "\\"));

        String outputName = "";
        if (type.equals("globalHeight")) {
            outputName = "globalHeight";
        } else if (type.equals("updateLocalHeight")) {
            outputName = "updateLocalHeight";
        } else if (type.equals("updateGlobalHeight")) {
            outputName = "updateGlobalHeight";
        }

        String finalOutputName1 = outputName;
        File[] files = dir.listFiles((dir1, name) -> (name.endsWith(".png") && name.startsWith(finalOutputName1)));

        BufferedImage first = ImageIO.read(new File(currentWorkingPath.concat("\\" + storageLocationName + "\\" + iterationLocation + "\\" + outputName + "_0" + ".png")));
        ImageOutputStream output = new FileImageOutputStream(new File(currentWorkingPath.concat("\\" + storageLocationName + "\\" + iterationLocation +
                "\\gif_" + outputName + ".gif")));

        GifSequenceWriter writer = new GifSequenceWriter(output, first.getType(), GIF_DELAY, true);
        writer.writeToSequence(first);


        String finalOutputName = outputName;
        Arrays.sort(files, new Comparator<File>() {
            @Override
            public int compare(File f1, File f2) {
                String s1 = f1.getName().substring(finalOutputName.length() + 1, f1.getName().indexOf("."));
                String s2 = f2.getName().substring(finalOutputName.length() + 1, f2.getName().indexOf("."));
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

    public static void applyHeightUpdate(Cell[][] grid, double[][] heightUpdate) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[i][j].height = grid[i][j].height + heightUpdate[i][j];

            }
        }

    }


    public static double[][] addHeights(double[][] previousHeight, double[][] newHeight) {

        double[][] result = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                // MEMORY_DECAY_RATE = 0.5; // ##hardcoded
                result[i][j] = previousHeight[i][j] * MEMORY_DECAY_RATE + newHeight[i][j];

            }
        }

        return result;
    }

    // adds the propagated height to the base and stores it in global grid
    public static Cell[][] addToBaseFunction(Cell[][] grid, double[][] base, double[][] previousHeight) {

        double[][] result = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[i][j].height = base[i][j] + previousHeight[i][j];

            }
        }

        return grid;
    }


    public static void sampleFlowAccumulation(Cell[][] grid) {

    }

    public static void computeFlowAccumulation(Cell[][] grid) {

        for (int col = 0; col < NR_OF_COLUMNS; col++) {

            for (int row = 0; row < NR_OF_ROWS; row++) {
                if (col == sourceX && row == sourceY) {
                    continue;
                }
                grid[col][row].flowAccumulation = 0;

            }
        }

        for (int col = 0; col < NR_OF_COLUMNS; col++) {

            for (int row = 0; row < NR_OF_ROWS; row++) {
                if (col == sourceX && row == sourceY) {
                    continue;
                }
                int flowDirection = grid[col][row].flowDirection;

                Cell flowsInto = null;

                if (flowDirection == 1) {
                    flowsInto = grid[col + 1][row];
                } else if (flowDirection == 2) {
                    flowsInto = grid[col + 1][row + 1];
                } else if (flowDirection == 4) {
                    flowsInto = grid[col][row + 1];
                } else if (flowDirection == 8) {
                    flowsInto = grid[col - 1][row + 1];
                } else if (flowDirection == 16) {
                    flowsInto = grid[col - 1][row];
                } else if (flowDirection == 32) {
                    flowsInto = grid[col - 1][row - 1];
                } else if (flowDirection == 64) {
                    flowsInto = grid[col][row - 1];
                } else if (flowDirection == 128) {
                    flowsInto = grid[col + 1][row - 1];
                }

                flowsInto.flowAccumulation = flowsInto.flowAccumulation + 1;

            }
        }

    }

    public static Tuple<Cell[][], ArrayList<Path>> iterate(
            Cell[][] grid, ArrayList pointsList, ArrayList paths,
            ArrayList distancesForPaths, int NR_OF_ITERATIONS,
            String BASE_HEIGHT_TYPE, double baseFunctionScale, boolean saveOutputs,
            boolean showIntermediateResults, double width, double scale, String iterationLocation,
            double[][] baseFunction)
            throws IOException {

        double[][] sumOfPreviousHeights = new double[NR_OF_ROWS][NR_OF_COLUMNS];

        if (MEMORY_MODE) {
            // this is the first height update w.r.t. the base

            // take the old height (that is, the global grid height), compute height update U_1
            double[][] heightUpdate = computeHeightUpdate(grid, distancesForPaths, paths, 0, iterationLocation); // this is U_1

            // here we apply the height update on the base function
            // the update is written inside global grid height
            applyHeightUpdate(grid, heightUpdate); // grid = B + U_1

            // sum of previous heights is now U_1
            sumOfPreviousHeights = heightUpdate;

        } else {
            updateHeightOverwrite(grid, distancesForPaths, width, scale, paths, 0, iterationLocation);
        }

        computeMinAndMaxHeights(grid);

        //adjustHeightMinDistances(grid, distancesForPaths, width, scale, paths);

        // Iterate a number of times:
        for (int iteration = 1; iteration < NR_OF_ITERATIONS; iteration++) {
            System.out.println("iteration : " + iteration);
            logFileWriter.write("iteration : " + iteration + "\n");

            computeFlow(grid, iteration);

            // computeFlowAccumulation(grid);
            ArrayList<Path> oldPaths = paths;

            if (HORIZONTAL_FLOW_MODE) {
                paths = computePathsToFrameEdge(pointsList, grid);
            } else {
                paths = computePaths(pointsList, grid);
            }
            boolean converged = testConvergence(oldPaths, paths);
            System.out.println("Converged : " + converged);

            if (paths == null) {
                return null;
            }

            if (saveOutputs == true) {
                if (GENERATE_INTERMEDIATE_RESULTS) {
                    draw(grid, paths, iteration, showIntermediateResults, iterationLocation, width, scale, pointsList);
                    //drawPaths(grid, paths, i, showIntermediateResults, iterationLocation, width, scale);
                    //drawFlow(grid, paths, i, showIntermediateResults, iterationLocation, width, scale);
                    //drawFlowAccumulation(grid, paths, i, false, iterationLocation, width, scale);

                }

                if (DRAW_DISTANCE_IMAGES) {
                    drawDistances(grid, paths, distancesForPaths, showIntermediateResults, width, scale, iteration, iterationLocation);
                }
            }
            if (converged) {
                //return null;
            }

            long startTime = System.currentTimeMillis();

            if (DISTANCE_METRIC.equals("BFS")) {
                distancesForPaths = computeBfs(grid, paths);
            } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                distancesForPaths = compute_Dijkstra(grid, paths);
            } else if (DISTANCE_METRIC.equals("ARC")) {
                distancesForPaths = computeArcLength(grid, paths);
            } else if (DISTANCE_METRIC.equals("ANGULAR_INTERSECTION")) {
                distancesForPaths = computeAngularDistanceWithIntersection(grid, paths);
            } else if (DISTANCE_METRIC.equals("ANGULAR_WITH_ARC_LENGTH")) {
                distancesForPaths = computeAnguarWithArcLength(grid, paths);
            } else if (DISTANCE_METRIC.equals("POLAR_SYSTEM")) {
                distancesForPaths = computeVerticalDistances(grid, paths);
            } else if (DISTANCE_METRIC.equals("BRUTE_FORCE")) {
                distancesForPaths = computeDistancesBruteForce(grid, paths);
            }

            long endTime = System.currentTimeMillis();

            Iterator pathIt = distancesForPaths.iterator();

            ArrayList transposedDist = new ArrayList();

            while (pathIt.hasNext()) {

                DistanceForPathMatrix dist = (DistanceForPathMatrix) pathIt.next();
                double[][] distMatrix = dist.distanceMatrix;

                transposeMatrix(distMatrix);
                DistanceForPathMatrix distTransposed = new DistanceForPathMatrix();
                distTransposed.distanceMatrix = distMatrix;
                distTransposed.pathId = dist.pathId;
                transposedDist.add(distTransposed);

            }

            distancesForPaths = transposedDist;
            System.out.println("That took " + (endTime - startTime) + " milliseconds");
            logFileWriter.write("That took " + (endTime - startTime) + " milliseconds" + "\n");

            if (RESET_HEIGHTS == true) {
                if (BASE_HEIGHT_TYPE.equals("EUCLID")) {
                    initializeGridHeight_EuclideanDist(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQUARED")) {
                    initializeGridHeight_EuclideanSquared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("chebyshev")) {
                    initializeGridHeightChebyshevDistance(grid);
                } else if (BASE_HEIGHT_TYPE.equals("EUCLID_SQRT")) {
                    initializeGridHeight_EuclideanDistSqrt(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE")) {
                    initializeGridHeightToEdge(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQUARED")) {
                    initializeGridHeightToEdgeSquared(grid);
                } else if (BASE_HEIGHT_TYPE.equals("TO_EDGE_SQRT")) {
                    initializeGridHeightToEdgeSqrt(grid);
                }
            }

            if (MEMORY_MODE) {
                // take the global grid height and compute the height update // here we took B + U_1 and computed U_2
                double[][] heightUpdate = computeHeightUpdate(grid, distancesForPaths, paths, iteration, iterationLocation); // = U_2

                sumOfPreviousHeights = addHeights(sumOfPreviousHeights, heightUpdate);

                // apply the height update on the global grid // here we should take B and add 1/2 U_1 + U_2
                //applyHeightUpdate(grid, heightUpdate); // this is not right

                grid = addToBaseFunction(grid, baseFunction, sumOfPreviousHeights);

            } else {
                updateHeightOverwrite(grid, distancesForPaths, width, scale, paths, iteration, iterationLocation);
            }

            double[][] heightUpdateObstacles;
            computeMinAndMaxHeights(grid);

            if (OBSTACLES.equals("STATIC")) {
                heightUpdateObstacles = computeObstaclesStatic(grid, pointsList);
                drawObstacles(heightUpdateObstacles, paths, iteration, iterationLocation);

            } else if (OBSTACLES.equals("PROGRESSIVE")) {
                heightUpdateObstacles = initializeObstaclesProgressive(grid, pointsList, iteration);
                drawObstacles(heightUpdateObstacles, paths, iteration, iterationLocation);

            }

           // computeMinAndMaxHeights(grid);

            //adjustHeightMinDistances(grid, distancesForPaths, width, scale, paths);

        }

        computeFlow(grid, NR_OF_ITERATIONS);
        // computeFlowAccumulation(grid);

        if (HORIZONTAL_FLOW_MODE) {
            paths = computePathsToFrameEdge(pointsList, grid);
        } else {
            paths = computePaths(pointsList, grid);
        }

        if (paths == null) {
            return null;
        }

        Tuple<Cell[][], ArrayList<Path>> tuple = new Tuple<Cell[][], ArrayList<Path>>(grid, paths);

        return tuple;

    }

    public static boolean testConvergence(ArrayList<Path> oldPaths, ArrayList<Path> newPaths) {

        for (int i = 0; i < oldPaths.size(); i++) {

            Path oldPath = oldPaths.get(i);
            Path newPath = newPaths.get(i);

            for (int j = 0; j < oldPath.cells.size(); j++) {

                Cell oldCell = (Cell) oldPath.cells.get(j);
                Cell newCell = (Cell) newPath.cells.get(j);

                if (oldCell.cellRow == newCell.cellRow && oldCell.cellCol == newCell.cellCol) {
                    // cells were the same
                } else {
                    System.out.println("did not converge");
                    return false;
                }

            }

        }
        System.out.println("Converged");
        return true;
    }

    public static ArrayList computeY_coordinatesForColumns(Cell[][] grid, ArrayList paths) {

        // TODO: potentially make this faster with binary search or something

        ArrayList yCoordinatesForAllColumns = new ArrayList();

        for (int col = 0; col < NR_OF_COLUMNS; col++) {

            ArrayList yCoordinatesForColumn = new ArrayList();

            for (int row = 0; row < NR_OF_ROWS; row++) {

                Iterator pathIterator = paths.iterator();
                int intersectionCounter = 0;

                while (pathIterator.hasNext()) {

                    Path path = (Path) pathIterator.next();

                    int pathId = path.id;

                    Iterator pathCellIterator = path.cells.iterator();

                    while (pathCellIterator.hasNext()) {

                        Cell pathCell = (Cell) pathCellIterator.next();

                        if (grid[col][row] == pathCell) {
                            intersectionCounter++;

                            // path cell is on the vertical cell
                            Pair<Integer, Integer> idAndY_coord = new Pair<>(pathId, row);

                            yCoordinatesForColumn.add(idAndY_coord);
                            //pathIndexes.add(pathIndex);
                        }
                    }
                }
            }

            //Map<Integer, List<Integer>> yCoordHashMap = new HashMap<>();

            // collects all the intersections and puts them in a list. Used to take the average of vertical intersection segments
            Map<Integer, List<Integer>> multiIntersections = new HashMap<Integer, List<Integer>>();
            for (int i = 0; i < yCoordinatesForColumn.size(); i++) {

                Pair pair = (Pair) yCoordinatesForColumn.get(i);

                List<Integer> list;
                if (multiIntersections.containsKey(pair.getKey())) {
                    list = multiIntersections.get(pair.getKey());
                    list.add((Integer) pair.getValue());
                } else {
                    list = new ArrayList<Integer>();
                    list.add((Integer) pair.getValue());
                    multiIntersections.put((Integer) pair.getKey(), list);
                }
            }

            for (Map.Entry<Integer, List<Integer>> set : multiIntersections.entrySet()) {

                int key = set.getKey();
                List<Integer> value = set.getValue();

                if (multiIntersections.get(key).size() > 1) {

                    int sum = 0;
                    for (int j = 0; j < multiIntersections.get(key).size(); j++) {
                        sum = sum + multiIntersections.get(key).get(j);
                    }
                    int avg = sum / multiIntersections.get(key).size();

                    if (multiIntersections.get(key).size() > 1) {
                        System.out.println();
                    }

                    Iterator yIterator = yCoordinatesForColumn.iterator();

                    ArrayList found = new ArrayList();
                    while (yIterator.hasNext()) {

                        Pair pair = (Pair) yIterator.next();

                        if ((Integer) pair.getKey() == key) {
                            found.add(pair);
                            //yCoordinatesForColumn.remove(pair);
                            //e--;
                        }

                    }
                    yCoordinatesForColumn.removeAll(found);

                    yCoordinatesForColumn.add(new Pair<>(key, avg));

                }

                // sort yCoordinatesForColumn
                //yCoordinatesForColumn.sort();

                // System.out.println();
                Collections.sort(yCoordinatesForColumn, new Comparator<Pair<Integer, Integer>>() {
                    @Override
                    public int compare(Pair<Integer, Integer> p1, Pair<Integer, Integer> p2) {
                        return p1.getValue().compareTo(p2.getValue());
                    }
                });
                // System.out.println();

            }


//            for (Map.Entry<Integer, List<Integer>> set : multiIntersections.entrySet()) {
//
////                try{
////                    if (multiIntersections.get(i).size() > 1);
////                } catch (Exception e) {
////                    System.out.println();
////                }
//                int i = set.getKey();
//                List<Integer> value = set.getValue();
//
//
//                if (multiIntersections.get(i).size() > 1) {
//
//                    int sum = 0;
//                    for (int j = 0; j < multiIntersections.get(i).size(); j++) {
//                        sum = sum + multiIntersections.get(i).get(j);
//                    }
//                    int avg = sum / multiIntersections.get(i).size();
//
//                    if (multiIntersections.get(i).size() > 3) {
//                        //System.out.println();
//                    }
//
//                    Iterator yIterator = yCoordinatesForColumn.iterator();
//
//                    ArrayList found = new ArrayList();
//                    while (yIterator.hasNext()) {
//
//                        Pair pair = (Pair) yIterator.next();
//
//                        if ((Integer) pair.getKey() == i) {
//                            found.add(pair);
//                            //yCoordinatesForColumn.remove(pair);
//                            //e--;
//                        }
//
//                    }
//                    yCoordinatesForColumn.removeAll(found);
//
//                    yCoordinatesForColumn.add(i, new Pair<>(i, avg));
//
//                }
//            }

            yCoordinatesForAllColumns.add(yCoordinatesForColumn);

        }

        return yCoordinatesForAllColumns;
    }

    public static double computeScaledPathFactors(Cell[][] grid, ArrayList paths, ArrayList distancesForCell, Cell cell) {

        ArrayList yCoordinates = new ArrayList();
        ArrayList pathIndexes = new ArrayList();

        for (int i = 0; i < NR_OF_ROWS; i++) {

            Iterator pathIterator = paths.iterator();

            int pathIndex = 0;

            while (pathIterator.hasNext()) {

                ArrayList path = (ArrayList) pathIterator.next();

                Iterator pathCellIterator = path.iterator();

                while (pathCellIterator.hasNext()) {

                    Cell pathCell = (Cell) pathCellIterator.next();

                    if (grid[cell.cellCol][i] == pathCell) {

                        // path cell is on the vertical cell
                        yCoordinates.add(i);
                        pathIndexes.add(pathIndex);
                    }
                }
            }
        }

        ArrayList factors = new ArrayList();

        int index = -1;

        for (int i = 0; i < yCoordinates.size(); i++) {

            if (i == yCoordinates.size() - 1) {
                index = i;
                break;
            }

            if (cell.cellRow <= (int) yCoordinates.get(i + 1)) {
                index = i;
                break;
            }

            if (cell.cellRow <= (int) yCoordinates.get(i + 1) && cell.cellRow >= (int) yCoordinates.get(i)) {
                index = i;
                break;
            }
        }

        for (int i = 0; i < yCoordinates.size(); i++) {

            int distanceFromIndex = Math.abs(i - index);
            factors.add(1 / (Math.pow(2, distanceFromIndex)));

        }

        double sum = 0;
        for (int k = 0; k < distancesForCell.size(); k++) {
            sum = sum + (double) factors.get(k) * gaussian((double) distancesForCell.get(k), 0, HEIGHT_FUNCTION_WIDTH);

        }

        return sum;
    }

    public static double gaussian(double x, double mu, double sigma) {

        return (double) (Math.exp(-Math.pow((x - mu) / sigma, 2.0) / 2));

    }

    public static double gaussian_2(double x, double mu, double sigma) {

        return (double) (Math.exp(-Math.pow((x - mu) / sigma, 2.0) / 2));

    }

    public static double heightFunction(double x, double width) {

        double result;

        double sigma = 0.5 * width; //width;//2 * width / 5;

        if (x < sigma) {
            result = Math.exp(((Math.pow(-x, 2)) / (2 * Math.pow(sigma, 2))));
        } else if (x < 2 * sigma) {
            result = Math.exp(-1 / 2) * (2 - x / sigma);
        } else {
            result = 0;
        }
        return result;
    }

    public static void adjustHeightMinDistances(Cell[][] grid, ArrayList distancesForPaths, double width, double scale, ArrayList paths) throws IOException {
        System.out.println("adjusting height");
        logFileWriter.write("adjusting height" + "\n");

        double[][] dist = (double[][]) distancesForPaths.get(0);

        double[][] minDistances = new double[NR_OF_ROWS][NR_OF_COLUMNS];

        for (int i = 0; i < NR_OF_ROWS; i++) {
            for (int j = 0; j < NR_OF_COLUMNS; j++) {

                Iterator distIter = distancesForPaths.iterator();

                double minDistance = (double) dist[i][j];

                while (distIter.hasNext()) {

                    double[][] distancesForPath = (double[][]) distIter.next();

                    if (distancesForPath[i][j] < minDistance) {

                        minDistance = distancesForPath[i][j];

                    }

                }
                minDistances[i][j] = minDistance;
            }
        }

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                Cell cell = grid[j][i];

                double newHeight = gaussian((double) minDistances[cell.cellRow][cell.cellCol], 0, HEIGHT_FUNCTION_WIDTH);

                double height = ((-HEIGHT_FUNCTION_SCALE * newHeight));
                grid[j][i].height = grid[j][i].height + ((height));

            }
        }
    }

    public static Map<Pair, List<Cell>> computePathOverlaps(ArrayList<Path> paths) {

        // holds the info about overlapping cells.
        // indexed by pair of two overlapping paths. Value of index is the overlapping cell
        Map<Pair, List<Cell>> overlaps = new HashMap<Pair, List<Cell>>();

        for (int i = 0; i < paths.size(); i++) {

            Path path_1 = paths.get(i);

            for (int j = i + 1; j < paths.size(); j++) {

                if (i == j) {
                    continue;
                }

                Path path_2 = paths.get(j);

                Iterator path_1Iterator = path_1.cells.iterator();

                while (path_1Iterator.hasNext()) {

                    Cell path_1Cell = (Cell) path_1Iterator.next();

                    Iterator path_2Iterator = path_2.cells.iterator();

                    while (path_2Iterator.hasNext()) {

                        Cell path_2Cell = (Cell) path_2Iterator.next();

                        // TODO: potentiall check if they are within some range
                        if (path_1Cell.cellRow == path_2Cell.cellRow &&
                                path_1Cell.cellCol == path_2Cell.cellCol) {

                            // we have an overlapping cell between path_1 and path_2

                            Pair<Integer, Integer> overlappingPathIds = new Pair<>(i, j);
                            //overlaps.put(overlappingPathIds, path_1Cell);

                            List<Cell> list;

                            if (overlaps.containsKey(overlappingPathIds)) {
                                list = overlaps.get(overlappingPathIds);
                                list.add(path_1Cell);
                            } else {
                                list = new ArrayList<Cell>();
                                list.add(path_1Cell);
                                overlaps.put((Pair) overlappingPathIds, list);
                            }
                        }
                    }
                }
            }
        }

        if (overlaps.size() > 0) {
            System.out.println();

            for (HashMap.Entry<Pair, List<Cell>> entry_1 : overlaps.entrySet()) {

                Pair twoOverlappingPathIds_1 = entry_1.getKey();
                int hash1 = System.identityHashCode(twoOverlappingPathIds_1);
                List<Cell> listOfOverlappingCells_1 = entry_1.getValue();

                for (HashMap.Entry<Pair, List<Cell>> entry_2 : overlaps.entrySet()) {

                    Pair twoOverlappingPathIds_2 = entry_2.getKey();
                    if (entry_1 == entry_2) {
                        continue;
                    }
                    if (hash1 > System.identityHashCode(twoOverlappingPathIds_2)) {
                        continue;
                    }

                    List<Cell> listOfOverlappingCells_2 = entry_2.getValue();

                    if (listOfOverlappingCells_1.containsAll(listOfOverlappingCells_2) && !listOfOverlappingCells_2.containsAll(listOfOverlappingCells_1)) {
                        System.out.println();

                        listOfOverlappingCells_1.removeAll(listOfOverlappingCells_2);

                    }
                }
            }
        }


        HashMap<List<Cell>, Integer> overlapCount = new HashMap<>();

        for (HashMap.Entry<Pair, List<Cell>> entry_1 : overlaps.entrySet()) {

            Pair twoOverlappingPathIds_1 = entry_1.getKey();
            int hash1 = System.identityHashCode(twoOverlappingPathIds_1);
            List<Cell> listOfOverlappingCells_1 = entry_1.getValue();

            for (HashMap.Entry<Pair, List<Cell>> entry_2 : overlaps.entrySet()) {

                Pair twoOverlappingPathIds_2 = entry_2.getKey();

                if (entry_1 == entry_2) {
                    continue;
                }

                if (hash1 > System.identityHashCode(twoOverlappingPathIds_2)) {
                    continue;
                }

                List<Cell> listOfOverlappingCells_2 = entry_2.getValue();

                if (listOfOverlappingCells_1.containsAll(listOfOverlappingCells_2) && listOfOverlappingCells_2.containsAll(listOfOverlappingCells_1)) {
                    System.out.println();
                    // we have two equal paths

                    //Path path = new Path();
                    //path.cells = (ArrayList) listOfOverlappingCells_1;

                    if (overlapCount.containsKey(listOfOverlappingCells_1)) {
                        int count = overlapCount.get(listOfOverlappingCells_1);
                        overlapCount.put(listOfOverlappingCells_1, count + 1);
                    } else if (!overlapCount.containsKey(listOfOverlappingCells_1)) {
                        overlapCount.put(listOfOverlappingCells_1, 2);
                    }

                }
            }
        }

        return overlaps;

    }

    public static ArrayList computeMaxDistancesForColumns(ArrayList yCoordinatesOfPathsForCell) {

        ArrayList maxDistForAllCols = new ArrayList();

        for (int col = 0; col < NR_OF_COLUMNS; col++) {

            HashMap<Integer, Integer> maxY_neighborDistancesForColumn = new HashMap();
            for (int k = 0; k < yCoordinatesOfPathsForCell.size(); k++) {

                int distanceToFurthestNeighbor;

                if (k == 0) {

                    if (yCoordinatesOfPathsForCell.size() == 1) {
                        // if there is only one intersection in the column :
                        // we should take the largest distance to the edge
                        Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(0);

                        distanceToFurthestNeighbor = Math.abs((int) currentY.getValue() - NR_OF_COLUMNS);

                        distanceToFurthestNeighbor = Math.max((int) currentY.getValue(),
                                distanceToFurthestNeighbor);

                        maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);
                    } else if (yCoordinatesOfPathsForCell.size() > 1) {
                        Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(0);
                        Pair<Integer, Integer> nextY = (Pair) yCoordinatesOfPathsForCell.get(1);

                        distanceToFurthestNeighbor = Math.abs((int) currentY.getValue() - (int) nextY.getValue());

                        maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);
                    }

                } else if (k > 0 && k < yCoordinatesOfPathsForCell.size() - 1) {

                    Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(k);
                    Pair<Integer, Integer> nextY = (Pair) yCoordinatesOfPathsForCell.get(k + 1);
                    Pair<Integer, Integer> previousY = (Pair) yCoordinatesOfPathsForCell.get(k - 1);

                    distanceToFurthestNeighbor = Math.max(Math.abs((int) currentY.getValue() - (int) nextY.getValue()),
                            Math.abs((int) currentY.getValue() - (int) previousY.getValue()));

                    maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);

                } else if (k > 0 && k == yCoordinatesOfPathsForCell.size() - 1) {

                    Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(k);
                    Pair<Integer, Integer> previousY = (Pair) yCoordinatesOfPathsForCell.get(k - 1);

                    distanceToFurthestNeighbor = Math.abs((int) currentY.getValue() - (int) previousY.getValue());

                    maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);

                }
            }

            maxDistForAllCols.add(maxY_neighborDistancesForColumn);

        }


        return maxDistForAllCols;

    }

    public static double[][] computeHeightUpdate(Cell[][] grid, ArrayList distancesForPaths, ArrayList paths,
                                                 int iteration, String iterationLocation) throws IOException {

        System.out.println("adjusting height");
        logFileWriter.write("adjusting height" + "\n");

        double[][] computedHeight = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        ArrayList yCoordinatesForColumns = null;
        ArrayList maxDistancesForColumns = null;
        Map<Pair, List<Cell>> overlaps = null;
        if (PATH_SCALING) {
            if (SCALING_MODE.equals("WIDTHS")) {
                yCoordinatesForColumns = computeY_coordinatesForColumns(grid, paths);

                // TODO: compute max distances for paths here (for all columns)
                //maxDistancesForColumns = computeMaxDistancesForColumns(yCoordinatesForColumns);

                //System.out.println();

            }

            if (SCALING_MODE.equals("OVERLAPS")) {
                overlaps = computePathOverlaps(paths);

                ArrayList overlappingPaths = new ArrayList();

                for (HashMap.Entry<Pair, List<Cell>> entry_1 : overlaps.entrySet()) {

                    Pair twoOverlappingPathIds = entry_1.getKey();
                    ArrayList<Cell> listOfOverlappingCells = (ArrayList<Cell>) entry_1.getValue();

                    Path path = new Path();
                    path.cells = listOfOverlappingCells;

                    overlappingPaths.add(path);

                }

                HashMap<Path, Integer> overlappingPathAndNrOfOverlaps = new HashMap<>();


                if (SCALING_MODE.equals("OVERLAPS")) {
                    if (DISTANCE_METRIC.equals("BFS")) {
                        distancesForPaths = computeBfs(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                        distancesForPaths = compute_Dijkstra(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("ARC")) {
                        distancesForPaths = computeArcLength(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("ANGULAR_INTERSECTION")) {
                        distancesForPaths = computeAngularDistanceWithIntersection(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("ANGULAR_WITH_ARC_LENGTH")) {
                        distancesForPaths = computeAnguarWithArcLength(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("POLAR_SYSTEM")) {
                        distancesForPaths = computeVerticalDistances(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("BRUTE_FORCE")) {
                        distancesForPaths = computeDistancesBruteForce(grid, paths);
                    }
                }
            }
        }

        for (int row = 0; row < NR_OF_ROWS; row++) { // i is the row id
            //System.out.println();
            for (int col = 0; col < NR_OF_COLUMNS; col++) { //j is the column id

                if (PATH_SCALING) {

                    //SCALING_MODE = "FACTORS";

                    if (SCALING_MODE.equals("FACTORS")) {

                        factorsForPaths(grid, col, row, paths,
                                yCoordinatesForColumns, computedHeight,
                                distancesForPaths);

                    } else if (SCALING_MODE.equals("WIDTHS")) {

                        widthsForPaths_2(grid, col, row, paths,
                                yCoordinatesForColumns, computedHeight,
                                distancesForPaths);

                    } else if (SCALING_MODE.equals("OVERLAPS")) {

                        heightUpdateOverlaps(grid, col, row, paths,
                                yCoordinatesForColumns, computedHeight,
                                distancesForPaths, overlaps);

                    }


                    // =================== PATH SCALING = FALSE
                    // ===================
                    // ===================
                    // ===================

                } else {

                    Cell cell = grid[col][row];

                    Iterator pathIterator = distancesForPaths.iterator();

                    ArrayList distancesForCell = new ArrayList();

                    while (pathIterator.hasNext()) {

                        DistanceForPathMatrix distances = (DistanceForPathMatrix) pathIterator.next();
                        double[][] distancesMatrix = distances.distanceMatrix;

                        distancesForCell.add(distancesMatrix[cell.cellRow][cell.cellCol]);

                    }

                    double sum = 0;
                    for (int k = 0; k < distancesForCell.size(); k++) {

                        double distance = (double) distancesForCell.get(k);

                        sum = sum + gaussian(distance, 0, HEIGHT_FUNCTION_WIDTH);
                    }
                    //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

                    long startTime = System.currentTimeMillis();

                    //double sum_2 = computeScaledPathFactors(grid, paths, distancesForCell, cell);

                    long endTime = System.currentTimeMillis();
                    //System.out.println("That took " + (endTime - startTime) + " milliseconds");

                    if (RESET_HEIGHTS == false) {

                        // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        //grid[col][row].height = grid[col][row].height + ((height));
                        computedHeight[col][row] = height;

                    } else {
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        computedHeight[col][row] = height;

                        //grid[col][row].height = grid[col][row].height + ((height));
                    }
                }
            }
        }


        if (GENERATE_INTERMEDIATE_RESULTS) {
            if (GENERATE_INTERMEDIATE_HEIGHT) {
                computeMinAndMaxHeights(grid);
                drawMatrix(computedHeight, paths, iteration, iterationLocation, false);
                drawMatrix(computedHeight, paths, iteration, iterationLocation, true);

            }
        }

        return computedHeight;
    }

    public static void heightUpdateOverlaps(Cell[][] grid, int col, int row, ArrayList paths,
                                            ArrayList yCoordinatesForColumns, double[][] computedHeight,
                                            ArrayList distancesForPaths, Map<Pair, List<Cell>> overlaps) {


        Cell cell = grid[col][row];

        Iterator pathIterator = distancesForPaths.iterator();

        ArrayList distancesForCell = new ArrayList();

        while (pathIterator.hasNext()) {

            DistanceForPathMatrix distances = (DistanceForPathMatrix) pathIterator.next();
            double[][] distancesMatrix = distances.distanceMatrix;

            distancesForCell.add(distancesMatrix[cell.cellRow][cell.cellCol]);

        }

        double sum = 0;
        for (int k = 0; k < distancesForCell.size(); k++) {
            sum = sum + gaussian((double) distancesForCell.get(k), 0, HEIGHT_FUNCTION_WIDTH);
        }
        //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

        long startTime = System.currentTimeMillis();

        //double sum_2 = computeScaledPathFactors(grid, paths, distancesForCell, cell);

        long endTime = System.currentTimeMillis();
        //System.out.println("That took " + (endTime - startTime) + " milliseconds");

        if (RESET_HEIGHTS == false) {

            // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            //grid[col][row].height = grid[col][row].height + ((height));
            computedHeight[col][row] = height;

        } else {
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            computedHeight[col][row] = height;

            //grid[col][row].height = grid[col][row].height + ((height));
        }


    }

    public static void widthsForPaths(Cell[][] grid, int col, int row, ArrayList paths,
                                      ArrayList yCoordinatesForColumns, double[][] computedHeight,
                                      ArrayList distancesForPaths) {

        Cell cell = grid[col][row];

        // return the column of yCoordinate intersections of paths for the cell
        ArrayList yCoordinatesOfPathsForCell = (ArrayList) yCoordinatesForColumns.get(cell.cellCol);

        // computes the absolute y distance from cell to all the intersections
        HashMap absY_distancesForCell = new HashMap();
        for (int k = 0; k < yCoordinatesOfPathsForCell.size(); k++) {

            Pair pair = (Pair) yCoordinatesOfPathsForCell.get(k);

            //Pair absD_pair = new Pair(pair.getKey(), Math.abs(row - (Integer) pair.getValue()));

            absY_distancesForCell.put(pair.getKey(), Math.abs(row - (Integer) pair.getValue()));

        }

        //System.out.println();

        Iterator pathIterator = distancesForPaths.iterator();
        ArrayList distancesForCell = new ArrayList();

        while (pathIterator.hasNext()) {

            DistanceForPathMatrix distanceForPathMatrix = (DistanceForPathMatrix) pathIterator.next();

            Pair distancesForCellAndId = new Pair(distanceForPathMatrix.pathId, distanceForPathMatrix.distanceMatrix[cell.cellRow][cell.cellCol]);

            distancesForCell.add(distancesForCellAndId);

        }

        double sum = 0;

        // loops over all the distances for all paths (even the ones outside column)
        for (int k = 0; k < distancesForCell.size(); k++) {

            Pair distancesForCellAndId = (Pair) distancesForCell.get(k);

            double distanceForPath = (double) distancesForCellAndId.getValue();

            int factor = 1;

            // Check if the path is in the intersections column. If not, the factor is 1
            if (absY_distancesForCell.containsKey(distancesForCellAndId.getKey())) {
                factor = ((int) absY_distancesForCell.get(distancesForCellAndId.getKey()) + 1);
            }

            sum = sum + gaussian((double) distanceForPath, 0, HEIGHT_FUNCTION_WIDTH);

        }

        if (RESET_HEIGHTS == false) {

            // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));   // grid[col][row]
            grid[col][row].height = grid[col][row].height + ((height)); // updates horizontally
            computedHeight[col][row] = height;

        } else {
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            computedHeight[col][row] = height;

            grid[col][row].height = grid[col][row].height + ((height));
        }

    }

    public static void widthsForPathsRadial() {


    }

    public static void widthsForPaths_2(Cell[][] grid, int col, int row, ArrayList paths,
                                        ArrayList yCoordinatesForColumns, double[][] computedHeight,
                                        ArrayList distancesForPaths) {

        Cell cell = grid[col][row];


        // here we compute the max distance to adjacent neighbors for each y intersection. This is indexed by the path id of the intersection.
        // TODO: this can be computed once per column! (currently it is computed for each cell)
//        HashMap<Integer, Integer> maxY_neighborDistancesForColumn = new HashMap();
//        for (int k = 0; k < yCoordinatesOfPathsForCell.size(); k++) {
//
//            int distanceToFurthestNeighbor;
//
//            if (k == 0) {
//
//                if (yCoordinatesOfPathsForCell.size() == 1) {
//                    // if there is only one intersection in the column :
//                    // we should take the largest distance to the edge
//                    Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(0);
//
//                    distanceToFurthestNeighbor = Math.abs((int) currentY.getValue() - NR_OF_COLUMNS);
//
//                    distanceToFurthestNeighbor = Math.max((int) currentY.getValue(),
//                            distanceToFurthestNeighbor);
//
//                    maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);
//                } else if (yCoordinatesOfPathsForCell.size() > 1) {
//                    Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(0);
//                    Pair<Integer, Integer> nextY = (Pair) yCoordinatesOfPathsForCell.get(1);
//
//                    distanceToFurthestNeighbor = Math.abs((int) currentY.getValue() - (int) nextY.getValue());
//
//                    maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);
//                }
//
//            } else if (k > 0 && k < yCoordinatesOfPathsForCell.size() - 1) {
//
//                Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(k);
//                Pair<Integer, Integer> nextY = (Pair) yCoordinatesOfPathsForCell.get(k + 1);
//                Pair<Integer, Integer> previousY = (Pair) yCoordinatesOfPathsForCell.get(k - 1);
//
//                distanceToFurthestNeighbor = Math.max(Math.abs((int) currentY.getValue() - (int) nextY.getValue()),
//                        Math.abs((int) currentY.getValue() - (int) previousY.getValue()));
//
//                maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);
//
//            } else if (k > 0 && k == yCoordinatesOfPathsForCell.size() - 1) {
//
//                Pair<Integer, Integer> currentY = (Pair) yCoordinatesOfPathsForCell.get(k);
//                Pair<Integer, Integer> previousY = (Pair) yCoordinatesOfPathsForCell.get(k - 1);
//
//                distanceToFurthestNeighbor = Math.abs((int) currentY.getValue() - (int) previousY.getValue());
//
//                maxY_neighborDistancesForColumn.put(currentY.getKey(), distanceToFurthestNeighbor);
//
//            }
//        }

        // return the column of yCoordinate intersections of paths for the cell
        ArrayList yCoordinatesOfPathsForCell = (ArrayList) yCoordinatesForColumns.get(cell.cellCol);

        //int width = 0;

        double update = 0.0;
        double constant = 0.5;//0.5;          //3
        double updateFactor = 1;
        HEIGHT_FUNCTION_SCALE = 1;     //3

        boolean pathOverlapped = false;

        //boolean circularMode = true;

        // loop over all paths
        for (int i = 0; i < yCoordinatesOfPathsForCell.size(); i++) {

            // intersection of path i
            Pair intersectionForPath_1 = (Pair) yCoordinatesOfPathsForCell.get(i);

            double width_1 = 0;
            if (i < yCoordinatesOfPathsForCell.size() - 1) {
                // get the intersection of next path if there is one. Width = distance from i to i + 1 (width of the rectangle defined by two intersections)
                width_1 = (int) ((Pair) yCoordinatesOfPathsForCell.get(i + 1)).getValue() - (int) intersectionForPath_1.getValue();
            } else {
                // if there is no next path, we loop to first path, width = total (looped) width of these two sections
                width_1 = NR_OF_ROWS - (int) intersectionForPath_1.getValue() + (int) ((Pair) yCoordinatesOfPathsForCell.get(0)).getValue();
            }

            double width_2 = 0.0;
            if (i > 0) {
                // if there is a path before path_1, we also compute width between i and i - 1
                width_2 = (int) intersectionForPath_1.getValue() - (int) ((Pair) yCoordinatesOfPathsForCell.get(i - 1)).getValue();
            } else {
                // if there was no path, we take the last path. Width = looped width of the two sections
                width_2 = NR_OF_ROWS + (int) intersectionForPath_1.getValue() -
                        (int) ((Pair) yCoordinatesOfPathsForCell.get(yCoordinatesOfPathsForCell.size() - 1)).getValue();
            }

            if (width_1 < 1 || width_2 < 1) {
                continue;
                //System.out.println();;
            }

            // width_1 and width_2 are widths of paths i - 1, i, i + 1. They are looped if out of bounds

            if (cell.cellRow == 32 && cell.cellCol == 0) {

                System.out.println();
            }
            if (cell.cellRow == 450 && cell.cellCol == 0) {
                System.out.println();

            }
            if (cell.cellRow == 470 && cell.cellCol == 0) {
                System.out.println();

            }

//            width_1 *= constant;
//            width_2 *= constant;

            // take the max width:
            //width_1 = Math.max(width_1, width_2);
            //width_1 = (width_1 + width_2) / 2;//Math.(width_1, width_2);
            //width_2 = width_1;

            // if path is above cell
            if (cell.cellRow > (int) intersectionForPath_1.getValue()) {

                // distance_1 = distance from current cell to intersection of path i
                double distance_1 = cell.cellRow - (int) intersectionForPath_1.getValue();
                // distance_2 = residual distance (everything except distance from cell to pI)
                double distance_2 = NR_OF_ROWS - distance_1;

                // influence_1 = contribution of path i
//                double influence_1 = gaussian_2(distance_1, 0, width_1);
//                // influence_2 = contribution of rest
//                double influence_2 = gaussian_2(distance_2, 0, width_2);
//

//                double normalizedDistance_1 = distance_1 / width_1;
//                double normalizedDistance_2 = distance_2 / width_2;

                double influence_1 = heightFunction(distance_1, width_1);
                double influence_2 = heightFunction(distance_2, width_2);

                update = update + ((influence_1 + influence_2)) * HEIGHT_FUNCTION_SCALE;

                //update = update + (influence) * HEIGHT_FUNCTION_SCALE * updateFactor;

            } else if (cell.cellRow < (int) intersectionForPath_1.getValue()) {
                // if path is below cell:
                // distance_1 = distance from path i above cell to cell
                double distance_1 = (int) intersectionForPath_1.getValue() - cell.cellRow;
                // distance_2 = residual distance
                double distance_2 = NR_OF_ROWS - distance_1;

//                double influence_1 = gaussian_2(distance_1, 0, width_2);
//                double influence_2 = gaussian_2(distance_2, 0, width_1);

//                double normalizedDistance_1 = distance_1 / width_1;
//                double normalizedDistance_2 = distance_2 / width_2;

                double influence_1 = heightFunction(distance_1, width_2);
                double influence_2 = heightFunction(distance_2, width_1);

                double influence = Math.max(influence_1, influence_2);

                update = update + ((influence_1 + influence_2)) * HEIGHT_FUNCTION_SCALE;
                //update = update + (influence) * HEIGHT_FUNCTION_SCALE * updateFactor;

            } else if (cell.cellRow == (int) intersectionForPath_1.getValue()) {

                // only add one overlap
                if (!pathOverlapped) {
                    update = update + (1 + heightFunction(NR_OF_ROWS, width_1)) * HEIGHT_FUNCTION_SCALE;
                    //update = update + (1 + gaussian_2(NR_OF_ROWS, 0, width_1)) * HEIGHT_FUNCTION_SCALE;

                    pathOverlapped = true;
                }
            }
        }

        if (RESET_HEIGHTS == false) {

            double height = ((-update));   // grid[col][row]
            grid[col][row].height = grid[col][row].height + ((height)); // updates horizontally
            computedHeight[col][row] = height;

        } else {
            double height = ((-update));
            computedHeight[col][row] = height;

            grid[col][row].height = grid[col][row].height + ((height));
        }

    }

    public static void factorsForPaths(Cell[][] grid, int col, int row, ArrayList paths,
                                       ArrayList yCoordinatesForColumns, double[][] computedHeight,
                                       ArrayList distancesForPaths) {

        Cell cell = grid[col][row];


        if (row == 5 && col == 5) {
            //  System.out.println();
        }

        if (row == 5 && col == 0) {
            // System.out.println();
        }

        if (row == 9 && col == 0) {
            //   System.out.println();
        }

        if (row == 5 && col == 8) {
            //    System.out.println();
        }

        // grid[col][row]
        // TODO: potentially has to be sorted s.t. the path ids are preserved
        ArrayList yCoordinatesOfPathsForCell = (ArrayList) yCoordinatesForColumns.get(cell.cellCol);

        //double[] factors = new double[paths.size()];

        Pair<Integer, Double>[] factorsForIds = new Pair[paths.size()];

        //int index = -1;

        // TODO: this could be binary search (if sorted)
        // Check where the cell is w.r.t. the intersections of paths (y coordinates)
        for (int p = 0; p < yCoordinatesOfPathsForCell.size(); p++) {

            Pair<Integer, Integer> pair = (Pair<Integer, Integer>) yCoordinatesOfPathsForCell.get(p);

            if (p + 1 < yCoordinatesOfPathsForCell.size()) {
                Pair<Integer, Integer> nextPair = (Pair<Integer, Integer>) yCoordinatesOfPathsForCell.get(p + 1);

                if (row > (int) pair.getValue() && row < (int) nextPair.getValue()) {
                    // index is between two y coordinates
                    // the two y coordinates should have 1

                    factorsForIds[p] = new Pair<>(pair.getKey(), 1.0);
                    factorsForIds[p + 1] = new Pair<>(nextPair.getKey(), 1.0);
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
                Pair<Integer, Integer> nextPair = (Pair<Integer, Integer>) yCoordinatesOfPathsForCell.get(p + counter);
                while (p + counter == nextPair.getValue()) {

                    factorsForIds[p + counter] = new Pair<>(nextPair.getKey(), 1.0);
                    //factors[p + counter] = 1;

                    counter++;

                    if (p + counter == yCoordinatesOfPathsForCell.size()) {
                        break;
                    }
                    nextPair = (Pair<Integer, Integer>) yCoordinatesOfPathsForCell.get(p + counter);
                }
                break;

            }

            if (p < pair.getValue()) {
                //index = p;
                factorsForIds[p] = new Pair<>(pair.getKey(), 1.0);
                //factors[p] = 1;
                break;
            }

        }

        int n = factorsForIds.length;
        int first = -1;
        int last = -1;

        for (int k = 0; k < n; k++) {
//                        if () {
//
//                        }

            if (factorsForIds[k] == null) {
                continue;
            }
            if (1.0 != factorsForIds[k].getValue()) {
                continue;
            }

            if (first == -1)
                first = k;
            last = k;
        }

        for (int p = 0; p < yCoordinatesOfPathsForCell.size(); p++) {

            int distanceFromIndexFirst = Math.abs(p - first);
            int distanceFromIndexLast = Math.abs(p - last);

            double factorFirst = 1 / (Math.pow(2, distanceFromIndexFirst));
            double factorLast = 1 / (Math.pow(2, distanceFromIndexLast));

            //Pair newPair = new Pair();
            Pair<Integer, Integer> pair = (Pair<Integer, Integer>) yCoordinatesOfPathsForCell.get(p);

            //try {
            factorsForIds[p] = new Pair<>(pair.getKey(), Math.max(factorFirst, factorLast));
//                        } catch (Exception e) {
//                            System.out.println();
//                        }

        }
        // System.out.println();

        Iterator pathIterator = distancesForPaths.iterator();
        ArrayList distancesForCell = new ArrayList();

        while (pathIterator.hasNext()) {

            DistanceForPathMatrix distanceForPathMatrix = (DistanceForPathMatrix) pathIterator.next();

            Pair distancesForCellAndId = new Pair(distanceForPathMatrix.pathId, distanceForPathMatrix.distanceMatrix[cell.cellRow][cell.cellCol]);

            distancesForCell.add(distancesForCellAndId);

        }

        double sum = 0;

        for (int k = 0; k < distancesForCell.size(); k++) {

            //System.out.println();
            Pair distancesForCellAndId = (Pair) distancesForCell.get(k);

            int pathId = (int) distancesForCellAndId.getKey();
            double distanceForPath = (double) distancesForCellAndId.getValue();
            double factor = 1;

            for (int b = 0; b < factorsForIds.length; b++) {

                if (factorsForIds[b] == null) {
                    factor = 1.0;
                    break;
                }

                if (factorsForIds[b].getKey() == k) {
                    factor = factorsForIds[b].getValue();
                    break;
                }

            }

            //double factorForPathId = factorsForIds[];

            //factor = 1.0;
            //System.out.print("row : " + row + " col : " + col + " factor : " + factor + " | ");
            //factor = 1.0;
            sum = sum + factor * gaussian((double) distanceForPath, 0, HEIGHT_FUNCTION_WIDTH);

//                        if (factors[k] == 0) {
//                            factors[k] = 1;
//                        }
//
//                        sum = sum + factors[k] * gaussian((double) distancesForCell.get(k), 0, HEIGHT_FUNCTION_WIDTH);
        }
        //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

        if (RESET_HEIGHTS == false) {

            // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));   // grid[col][row]
            grid[col][row].height = grid[col][row].height + ((height)); // updates horizontally
            computedHeight[col][row] = height;

        } else {
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            computedHeight[col][row] = height;

            grid[col][row].height = grid[col][row].height + ((height));
        }

    }

    public static void updateHeightOverwrite(Cell[][] grid, ArrayList distancesForPaths, double width, double scale, ArrayList paths, int iteration, String iterationLocation) throws IOException {

        System.out.println("adjusting height");
        logFileWriter.write("adjusting height" + "\n");

        boolean verbose = false;

        double[][] computedHeight = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        double[][] heightMatrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];

        for (int a = 0; a < NR_OF_COLUMNS; a++) {
            for (int b = 0; b < NR_OF_ROWS; b++) {

                heightMatrix[a][b] = grid[a][b].height;

            }
        }

        ArrayList yCoordinatesForColumns = null;
        ArrayList maxDistancesForColumns = null;
        Map<Pair, List<Cell>> overlaps = null;
        if (PATH_SCALING) {
            if (SCALING_MODE.equals("WIDTHS")) {
                yCoordinatesForColumns = computeY_coordinatesForColumns(grid, paths);

                // TODO: compute max distances for paths here (for all columns)
                //maxDistancesForColumns = computeMaxDistancesForColumns(yCoordinatesForColumns);

                //System.out.println();

            }

            if (SCALING_MODE.equals("OVERLAPS")) {
                overlaps = computePathOverlaps(paths);

                ArrayList overlappingPaths = new ArrayList();

                for (HashMap.Entry<Pair, List<Cell>> entry_1 : overlaps.entrySet()) {

                    Pair twoOverlappingPathIds = entry_1.getKey();
                    ArrayList<Cell> listOfOverlappingCells = (ArrayList<Cell>) entry_1.getValue();

                    Path path = new Path();
                    path.cells = listOfOverlappingCells;

                    overlappingPaths.add(path);

                }

                HashMap<Path, Integer> overlappingPathAndNrOfOverlaps = new HashMap<>();


                if (SCALING_MODE.equals("OVERLAPS")) {
                    if (DISTANCE_METRIC.equals("BFS")) {
                        distancesForPaths = computeBfs(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("DIJKSTRA")) {
                        distancesForPaths = compute_Dijkstra(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("ARC")) {
                        distancesForPaths = computeArcLength(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("ANGULAR_INTERSECTION")) {
                        distancesForPaths = computeAngularDistanceWithIntersection(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("ANGULAR_WITH_ARC_LENGTH")) {
                        distancesForPaths = computeAnguarWithArcLength(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("POLAR_SYSTEM")) {
                        distancesForPaths = computeVerticalDistances(grid, overlappingPaths);
                    } else if (DISTANCE_METRIC.equals("BRUTE_FORCE")) {
                        distancesForPaths = computeDistancesBruteForce(grid, paths);
                    }
                }
            }
        }

        for (int row = 0; row < NR_OF_ROWS; row++) { // i is the row id
            //System.out.println();
            for (int col = 0; col < NR_OF_COLUMNS; col++) { //j is the column id

                if (PATH_SCALING) {

                    if (SCALING_MODE.equals("FACTORS")) {

                        factorsForPaths(grid, col, row, paths,
                                yCoordinatesForColumns, computedHeight,
                                distancesForPaths);

                    } else if (SCALING_MODE.equals("WIDTHS")) {

                        widthsForPaths_2(grid, col, row, paths,
                                yCoordinatesForColumns, computedHeight,
                                distancesForPaths);

                    }

                    // =================== PATH SCALING = FALSE
                    // ===================
                    // ===================
                    // ===================

                } else {
                    Cell cell = grid[col][row];

                    Iterator pathIterator = distancesForPaths.iterator();

                    ArrayList distancesForCell = new ArrayList();

                    while (pathIterator.hasNext()) {

                        DistanceForPathMatrix distances = (DistanceForPathMatrix) pathIterator.next();
                        double[][] distancesMatrix = distances.distanceMatrix;

                        distancesForCell.add(distancesMatrix[cell.cellRow][cell.cellCol]);

                    }
                    double sum = 0;

                    if (FLOW_ACCUMULATION) {
//                        int flowAccumulationForCell = grid[col][row].flowAccumulation;
//                        double flowAccNormalized = (flowAccumulationForCell - minFlowAccumulation) / (minFlowAccumulation - maxFlowAccumulation);
//                        double factor = 1 - flowAccNormalized;
//                        factor = 1;
//
//                        for (int k = 0; k < distancesForCell.size(); k++) {
//
//                            double distance = (double) distancesForCell.get(k);
//
//                            sum = sum + gaussian(distance, 0, HEIGHT_FUNCTION_WIDTH) * factor;
//                        }
                    } else {
                        for (int k = 0; k < distancesForCell.size(); k++) {

                            double distance = (double) distancesForCell.get(k);
                            //sum = sum + heightFunction(distance, HEIGHT_FUNCTION_WIDTH);
                            sum = sum + gaussian(distance, 0, HEIGHT_FUNCTION_WIDTH);
                        }
                    }


                    //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

                    long startTime = System.currentTimeMillis();

                    //double sum_2 = computeScaledPathFactors(grid, paths, distancesForCell, cell);

                    long endTime = System.currentTimeMillis();
                    //System.out.println("That took " + (endTime - startTime) + " milliseconds");

                    if (RESET_HEIGHTS == false) {

                        // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        grid[col][row].height = grid[col][row].height + ((height));
                        computedHeight[col][row] = height;

                    } else {
                        double height = ((-HEIGHT_FUNCTION_SCALE * sum));
                        computedHeight[col][row] = height;
                        HEIGHT_UPDATE[col][row] = height;

                        grid[col][row].height = grid[col][row].height + ((height));
                    }
                }
            }
        }

        if (GENERATE_INTERMEDIATE_RESULTS) {
            if (GENERATE_INTERMEDIATE_HEIGHT) {
                computeMinAndMaxHeights(grid);
                drawMatrix(computedHeight, paths, iteration, iterationLocation, false);
                drawMatrix(computedHeight, paths, iteration, iterationLocation, true);

            }
        }
    }

    private static void heightUpdateOverlaps_2(Cell[][] grid, int col, int row, ArrayList paths,
                                               ArrayList yCoordinatesForColumns, double[][] computedHeight,
                                               ArrayList distancesForPaths, Map<Pair, List<Cell>> overlaps) {

        Cell cell = grid[col][row];


        for (int i = 0; i < paths.size(); i++) {


        }

        Iterator pathIterator = distancesForPaths.iterator();

        ArrayList distancesForCell = new ArrayList();

        while (pathIterator.hasNext()) {

            DistanceForPathMatrix distances = (DistanceForPathMatrix) pathIterator.next();
            double[][] distancesMatrix = distances.distanceMatrix;

            distancesForCell.add(distancesMatrix[cell.cellRow][cell.cellCol]);

        }

        double sum = 0;
        for (int k = 0; k < distancesForCell.size(); k++) {
            sum = sum + gaussian((double) distancesForCell.get(k), 0, HEIGHT_FUNCTION_WIDTH);
        }
        //System.out.print("i : " + row + " j : " + col + " : " + sum + " | ");

        long startTime = System.currentTimeMillis();

        //double sum_2 = computeScaledPathFactors(grid, paths, distancesForCell, cell);

        long endTime = System.currentTimeMillis();
        //System.out.println("That took " + (endTime - startTime) + " milliseconds");

        if (RESET_HEIGHTS == false) {

            // wtf is this??
//                    if (i == 0) {
//                        i = 1;
//                    }
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            //grid[col][row].height = grid[col][row].height + ((height));
            computedHeight[col][row] = height;

        } else {
            double height = ((-HEIGHT_FUNCTION_SCALE * sum));
            computedHeight[col][row] = height;

            //grid[col][row].height = grid[col][row].height + ((height));
        }

    }

    public static void initializeGridHeight_EuclideanSquared(Cell[][] grid) {

        double min = Integer.MAX_VALUE;
        double max = Integer.MIN_VALUE;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double height =  (Math.pow(grid[j][i].cellCol - sourceX, 2) + Math.pow(grid[j][i].cellRow - sourceY, 2));
                //height = (height - min) / (max - min);
                //grid[j][i].height = height;

                if (height < min) {
                    min = height;
                }
                if (height > max) {
                    max = height;
                }
            }
        }
        System.out.println("min: " + min);
        System.out.println("max: " + max);
        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double height =  (Math.pow(grid[j][i].cellCol - sourceX, 2) + Math.pow(grid[j][i].cellRow - sourceY, 2));
                height = 100 * (height - min) / (max - min);
                grid[j][i].height = height;
            }
        }

    }

    public static void initializeGridHeight_EuclideanDistSqrt(Cell[][] grid) {


        double min = Integer.MAX_VALUE;
        double max = Integer.MIN_VALUE;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double height = Math.sqrt(Math.sqrt(Math.pow(grid[j][i].cellCol - sourceX, 2) + Math.pow(grid[j][i].cellRow - sourceY, 2)));
                //height = (height - min) / (max - min);
                //grid[j][i].height = height;

                if (height < min) {
                    min = height;
                }
                if (height > max) {
                    max = height;
                }
            }
        }
        System.out.println("min: " + min);
        System.out.println("max: " + max);

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double height = Math.sqrt(Math.sqrt(Math.pow(grid[j][i].cellCol - sourceX, 2) + Math.pow(grid[j][i].cellRow - sourceY, 2)));

                height = 100 * (height - min) / (max - min);

                grid[j][i].height = height;

            }
        }

    }

    public static void initializeGridHeight_EuclideanDistSquared(Cell[][] grid) {

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = (BASE_SCALE * Math.pow(2, Math.sqrt(Math.pow(grid[j][i].cellCol - sourceX, 2) + Math.pow(grid[j][i].cellRow - sourceY, 2))));

            }
        }

    }

    public static void initializeGridHeightChebyshevDistance(Cell[][] grid) {

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


        Cell sourceCell = grid[sourceY][sourceX];

        // add all cells of a path to queue
        queue.add(sourceCell);

        visited[sourceCell.cellRow][sourceCell.cellCol] = true;

        while (!queue.isEmpty()) {

            Cell cell = queue.peek();

            int x = cell.cellCol;
            int y = cell.cellRow;

            queue.remove();

            for (int i = 0; i < 4; i++) {

                int adjX = x + dCol[i];
                int adjY = y + dRow[i];

                if (isValid(visited, adjY, adjX)) {

                    queue.add(grid[adjX][adjY]);
                    visited[adjY][adjX] = true;

                    distances[adjY][adjX] = distances[y][x] + 1;

                }
            }
        }

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                grid[j][i].height = distances[j][i];

            }
        }
    }


    public static void initializeGridHeight_EuclideanDist(Cell[][] grid) {

        double min = Integer.MAX_VALUE;
        double max = Integer.MIN_VALUE;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double height = (Math.sqrt(Math.pow(grid[j][i].cellCol - sourceX, 2) + Math.pow(grid[j][i].cellRow - sourceY, 2)));
                //height = (height - min) / (max - min);
                //grid[j][i].height = height;

                if (height < min) {
                    min = height;
                }
                if (height > max) {
                    max = height;
                }
            }
        }
        System.out.println("min before normalization: " + min);
        System.out.println("max before normalization: " + max);

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double height = (Math.sqrt(Math.pow(grid[j][i].cellCol - sourceX, 2) + Math.pow(grid[j][i].cellRow - sourceY, 2)));

                height = 100 * (height - min) / (max - min);

                grid[j][i].height = height;
            }
        }

        min = Integer.MAX_VALUE;
        max = Integer.MIN_VALUE;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double height = grid[i][j].height;
                //height = (height - min) / (max - min);
                //grid[j][i].height = height;

                if (height < min) {
                    min = height;
                }
                if (height > max) {
                    max = height;
                }
            }
        }
        System.out.println("min after normalization: " + min);
        System.out.println("max after normalization: " + max);

    }

    public static void assignZerosToGrid(Cell[][] grid) {

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

    public static boolean isInsideGrid(int x, int y) {

        if (x >= 0 && x < NR_OF_COLUMNS && y >= 0 && y < NR_OF_ROWS) {
            return true;
        } else {
            return false;
        }
    }

    public static ArrayList computeArcLength(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing arc length distance");
        logFileWriter.write("computing arc length distance" + "\n");

        Iterator pathIterator = paths.iterator();

        ArrayList distancesForPaths = new ArrayList();

        // for all paths
        while (pathIterator.hasNext()) {

            ArrayList path = (ArrayList) pathIterator.next();

            Collections.reverse(path);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(sourceCell.cellCol - cell.cellCol, 2) +
                            Math.pow(sourceCell.cellRow - cell.cellRow, 2)));

                    if (radius > ARC_RADIUS) {

                        // Binary search finds the first cell for which the distance to this cell is >= than radius

//                        if (i == 389 && j == 0) {
//                            System.out.println();
//                        }

                        int indexOfCell = binarySearch_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersectionCell = null;
                        //try {
                        intersectionCell = (Cell) path.get(indexOfCell);

//                        } catch (Exception e) {
//                            System.out.println();
//                        }

                        double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - intersectionCell.cellRow, 2)));

                        // TODO: check index out of bounds
                        if (indexOfCell + 1 < path.size() && indexOfCell - 1 > 0) {

                            Cell nextCell = (Cell) path.get(indexOfCell + 1);
                            Cell previousCell = (Cell) path.get(indexOfCell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(sourceCell.cellCol - nextCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - nextCell.cellRow, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(sourceCell.cellCol - previousCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - previousCell.cellRow, 2)));

                            Tuple<Double, Double> intersectionPoint = null;
                            // if radius is between intersectionCell and intersectionCell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        previousCell.cellCol, previousCell.cellRow);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersectionCell and intersectionCell + 1 we consider these two cells

                                // here compute the intersection

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        nextCell.cellCol, nextCell.cellRow);

                            } else if (radius == dist) {
                                intersectionPoint = new Tuple<Double, Double>((double) intersectionCell.cellCol, (double) intersectionCell.cellRow);
                            }

                            if (intersectionPoint == null) {
                                //System.out.println();
                            }

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(cell.cellRow - intersectionPoint.second, 2)));

                            double distanceFromSourceToIntersection = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(sourceCell.cellRow - intersectionPoint.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distanceFromSourceToIntersection, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * distanceFromSourceToIntersection));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength;

                        } else {

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionCell.cellCol, 2) +
                                    Math.pow(cell.cellRow - intersectionCell.cellRow, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * radius));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength;

                        }

                    } else {

                        int indexOfCell = binarySearch_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersectionCell = (Cell) path.get(indexOfCell);

                        double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - intersectionCell.cellRow, 2)));

                        // TODO: check index out of bounds
                        if (indexOfCell + 1 < path.size() && indexOfCell - 1 > 0) {

                            Cell nextCell = (Cell) path.get(indexOfCell + 1);
                            Cell previousCell = (Cell) path.get(indexOfCell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(sourceCell.cellCol - nextCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - nextCell.cellRow, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(sourceCell.cellCol - previousCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - previousCell.cellRow, 2)));

                            Tuple<Double, Double> intersectionPoint = null;
                            // if radius is between intersectionCell and intersectionCell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        previousCell.cellCol, previousCell.cellRow);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersectionCell and intersectionCell + 1 we consider these two cells

                                // here compute the intersection

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        nextCell.cellCol, nextCell.cellRow);

                            } else if (radius == dist) {
                                intersectionPoint = new Tuple<Double, Double>((double) intersectionCell.cellCol, (double) intersectionCell.cellRow);
                            }

                            if (intersectionPoint == null) {
                                //System.out.println();
                            }

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(cell.cellRow - intersectionPoint.second, 2)));

                            double distanceFromSourceToIntersection = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(sourceCell.cellRow - intersectionPoint.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distanceFromSourceToIntersection, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * distanceFromSourceToIntersection));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength;

                        } else {

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionCell.cellCol, 2) +
                                    Math.pow(cell.cellRow - intersectionCell.cellRow, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * radius));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength;
                        }
                    }
                }
            }

            Iterator pathCellIter = path.iterator();

            while (pathCellIter.hasNext()) {

                Cell cell = (Cell) pathCellIter.next();

                distances[cell.cellCol][cell.cellRow] = 0.0;

            }

            double[][] transposedMatrix = transposeMatrix(distances);

            distancesForPaths.add(transposedMatrix);

        }
        return distancesForPaths;
    }

    public static ArrayList computeAngularWith_Dijkstra(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular distance with arc length ");
        logFileWriter.write("computing angular distance with arc length " + "\n");

        Iterator pathIterator = paths.iterator();

        ArrayList distancesForPaths = new ArrayList();

        // for all paths
        while (pathIterator.hasNext()) {

            ArrayList path = (ArrayList) pathIterator.next();

            Collections.reverse(path);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(sourceCell.cellCol - cell.cellCol, 2) +
                            Math.pow(sourceCell.cellRow - cell.cellRow, 2)));

                    if (radius > ARC_RADIUS) {

                        // Binary search finds the first cell for which the distance to this cell is >= than radius
                        int indexOfCell = binarySearch_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersectionCell = (Cell) path.get(indexOfCell);

                        double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - intersectionCell.cellRow, 2)));

                        // TODO: check index out of bounds
                        if (indexOfCell + 1 < path.size() && indexOfCell - 1 > 0) {

                            Cell nextCell = (Cell) path.get(indexOfCell + 1);
                            Cell previousCell = (Cell) path.get(indexOfCell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(sourceCell.cellCol - nextCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - nextCell.cellRow, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(sourceCell.cellCol - previousCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - previousCell.cellRow, 2)));

                            Tuple<Double, Double> intersectionPoint = null;
                            // if radius is between intersectionCell and intersectionCell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        previousCell.cellCol, previousCell.cellRow);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersectionCell and intersectionCell + 1 we consider these two cells

                                // here compute the intersection

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        nextCell.cellCol, nextCell.cellRow);

                            } else if (radius == dist) {
                                intersectionPoint = new Tuple<Double, Double>((double) intersectionCell.cellCol, (double) intersectionCell.cellRow);
                            }

                            if (intersectionPoint == null) {
                                //System.out.println();
                            }

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(cell.cellRow - intersectionPoint.second, 2)));

                            double distanceFromSourceToIntersection = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(sourceCell.cellRow - intersectionPoint.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distanceFromSourceToIntersection, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * distanceFromSourceToIntersection));

                            distances[i][j] = angle;

                        } else {

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionCell.cellCol, 2) +
                                    Math.pow(cell.cellRow - intersectionCell.cellRow, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * radius));

                            distances[i][j] = angle;

                        }

                    } else {

                        int indexOfCell = binarySearch_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersectionCell = (Cell) path.get(indexOfCell);

                        double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - intersectionCell.cellRow, 2)));

                        // TODO: check index out of bounds
                        if (indexOfCell + 1 < path.size() && indexOfCell - 1 > 0) {

                            Cell nextCell = (Cell) path.get(indexOfCell + 1);
                            Cell previousCell = (Cell) path.get(indexOfCell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(sourceCell.cellCol - nextCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - nextCell.cellRow, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(sourceCell.cellCol - previousCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - previousCell.cellRow, 2)));

                            Tuple<Double, Double> intersectionPoint = null;
                            // if radius is between intersectionCell and intersectionCell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        previousCell.cellCol, previousCell.cellRow);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersectionCell and intersectionCell + 1 we consider these two cells

                                // here compute the intersection

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        nextCell.cellCol, nextCell.cellRow);

                            } else if (radius == dist) {
                                intersectionPoint = new Tuple<Double, Double>((double) intersectionCell.cellCol, (double) intersectionCell.cellRow);
                            }

                            if (intersectionPoint == null) {
                                // System.out.println();
                            }

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(cell.cellRow - intersectionPoint.second, 2)));

                            double distanceFromSourceToIntersection = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(sourceCell.cellRow - intersectionPoint.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distanceFromSourceToIntersection, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * distanceFromSourceToIntersection));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength / radius;

                        } else {

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionCell.cellCol, 2) +
                                    Math.pow(cell.cellRow - intersectionCell.cellRow, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * radius));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength / radius;
                        }

                    }
                }
            }

            Iterator pathCellIter = path.iterator();

            while (pathCellIter.hasNext()) {

                Cell cell = (Cell) pathCellIter.next();

                distances[cell.cellCol][cell.cellRow] = 0.0;

            }

            distancesForPaths.add(transposeMatrix(distances));
        }
        return distancesForPaths;

    }

    public static ArrayList computeAnguarWithArcLength(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular distance with arc length ");
        logFileWriter.write("computing angular distance with arc length " + "\n");

        Iterator pathIterator = paths.iterator();

        ArrayList distancesForPaths = new ArrayList();

        // for all paths
        while (pathIterator.hasNext()) {

            Path path = (Path) pathIterator.next();

            ArrayList pathCells = path.cells;

            Collections.reverse(pathCells);

            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distanceMatrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];
            distanceForPathMatrix.pathId = path.id;
            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(sourceCell.cellCol - cell.cellCol, 2) +
                            Math.pow(sourceCell.cellRow - cell.cellRow, 2)));

                    if (false){//(radius > ARC_RADIUS) {

                        // Binary search finds the first cell for which the distance to this cell is >= than radius
                        int indexOfCell = binarySearch_2(pathCells, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersectionCell = (Cell) pathCells.get(indexOfCell);

                        double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - intersectionCell.cellRow, 2)));

                        // TODO: check index out of bounds
                        if (indexOfCell + 1 < pathCells.size() && indexOfCell - 1 > 0) {

                            Cell nextCell = (Cell) pathCells.get(indexOfCell + 1);
                            Cell previousCell = (Cell) pathCells.get(indexOfCell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(sourceCell.cellCol - nextCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - nextCell.cellRow, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(sourceCell.cellCol - previousCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - previousCell.cellRow, 2)));

                            Tuple<Double, Double> intersectionPoint = null;
                            // if radius is between intersectionCell and intersectionCell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        previousCell.cellCol, previousCell.cellRow);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersectionCell and intersectionCell + 1 we consider these two cells

                                // here compute the intersection

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        nextCell.cellCol, nextCell.cellRow);

                            } else if (radius == dist) {
                                intersectionPoint = new Tuple<Double, Double>((double) intersectionCell.cellCol, (double) intersectionCell.cellRow);
                            }

                            if (intersectionPoint == null) {
                                // System.out.println();
                            }

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(cell.cellRow - intersectionPoint.second, 2)));

                            double distanceFromSourceToIntersection = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(sourceCell.cellRow - intersectionPoint.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distanceFromSourceToIntersection, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * distanceFromSourceToIntersection));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = angle;

                        } else {

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionCell.cellCol, 2) +
                                    Math.pow(cell.cellRow - intersectionCell.cellRow, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * radius));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = angle;

                        }

                    } else {

                        int indexOfCell = binarySearch_2(pathCells, radius);//binarySearch(path, 0, path.size(), radius);

                        Cell intersectionCell = (Cell) pathCells.get(indexOfCell);

                        double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - intersectionCell.cellRow, 2)));

                        // TODO: check index out of bounds
                        if (indexOfCell + 1 < pathCells.size() && indexOfCell - 1 > 0) {

                            Cell nextCell = (Cell) pathCells.get(indexOfCell + 1);
                            Cell previousCell = (Cell) pathCells.get(indexOfCell - 1);

                            double dist_1 = (Math.sqrt(Math.pow(sourceCell.cellCol - nextCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - nextCell.cellRow, 2)));

                            double dist_2 = (Math.sqrt(Math.pow(sourceCell.cellCol - previousCell.cellCol, 2) +
                                    Math.pow(sourceCell.cellRow - previousCell.cellRow, 2)));

                            Tuple<Double, Double> intersectionPoint = null;
                            // if radius is between intersectionCell and intersectionCell - 1, we consider these two cells
                            if (radius > dist_2 && radius < dist) {

                                // here compute the intersection point

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        previousCell.cellCol, previousCell.cellRow);


                            } else if (radius > dist && radius < dist_1) {
                                // if radius is between intersectionCell and intersectionCell + 1 we consider these two cells

                                // here compute the intersection

                                intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                        sourceX, sourceY, radius,
                                        intersectionCell.cellCol, intersectionCell.cellRow,
                                        nextCell.cellCol, nextCell.cellRow);

                            } else if (radius == dist) {
                                intersectionPoint = new Tuple<Double, Double>((double) intersectionCell.cellCol, (double) intersectionCell.cellRow);
                            }

                            if (intersectionPoint == null) {
                                //System.out.println();
                            }

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(cell.cellRow - intersectionPoint.second, 2)));

                            double distanceFromSourceToIntersection = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionPoint.first, 2) +
                                    Math.pow(sourceCell.cellRow - intersectionPoint.second, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distanceFromSourceToIntersection, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * distanceFromSourceToIntersection));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength / ARC_RADIUS;

                        } else {

                            double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionCell.cellCol, 2) +
                                    Math.pow(cell.cellRow - intersectionCell.cellRow, 2)));

                            double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                    (2.0 * radius * radius));

                            double arcLength = 2 * Math.PI * radius * (angle / 360);

                            distances[i][j] = arcLength / ARC_RADIUS;
                        }
                    }
                }
            }

            Iterator pathCellIter = pathCells.iterator();

            while (pathCellIter.hasNext()) {

                Cell cell = (Cell) pathCellIter.next();

                distances[cell.cellCol][cell.cellRow] = 0.0;

            }

            double[][] transposedMatrix = transposeMatrix(distances);

            distanceForPathMatrix.distanceMatrix = transposedMatrix;
            distancesForPaths.add(distanceForPathMatrix);

//            for (int i = 0; i < NR_OF_COLUMNS; i++) {
//                for (int j = 0; j < NR_OF_ROWS; j++) {
//
//                    Cell cell = grid[i][j];
//
//                    double radius = (Math.sqrt(Math.pow(sourceCell.cellX - cell.cellX, 2) +
//                            Math.pow(sourceCell.cellY - cell.cellY, 2)));
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
        return distancesForPaths;

    }

    public static ArrayList computeAngularDistanceWithIntersection(Cell[][] grid, ArrayList<Path> paths) throws IOException {

        System.out.println("computing angular intersections distance");
        logFileWriter.write("computing angular intersections distance" + "\n");

        Iterator pathIterator = paths.iterator();

        ArrayList distancesForPaths = new ArrayList();


        // for all paths
        while (pathIterator.hasNext()) {

            Path path = (Path) pathIterator.next();

            ArrayList pathCells = path.cells;

            Collections.reverse(pathCells);

            //double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distanceMatrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];
            distanceForPathMatrix.pathId = path.id;

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(sourceCell.cellCol - cell.cellCol, 2) +
                            Math.pow(sourceCell.cellRow - cell.cellRow, 2)));

//                    if (radius <= 2) {
//                        distanceForPathMatrix.distanceMatrix[i][j] = 0.0;
//                        continue;
//                    }

                    int indexOfCell = binarySearch_2(pathCells, radius);//binarySearch(path, 0, path.size(), radius);

                    Cell intersectionCell = (Cell) pathCells.get(indexOfCell);

                    double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionCell.cellCol, 2) +
                            Math.pow(sourceCell.cellRow - intersectionCell.cellRow, 2)));

                    // TODO: check index out of bounds
                    if (indexOfCell + 1 < pathCells.size() && indexOfCell - 1 > 0) {
                        Cell nextCell = null;
                        //try {
                        nextCell = (Cell) pathCells.get(indexOfCell + 1);

                        //} catch (Exception e) {
                        //    System.out.println();
                        // }
                        Cell previousCell = (Cell) pathCells.get(indexOfCell - 1);

                        double dist_1 = (Math.sqrt(Math.pow(sourceCell.cellCol - nextCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - nextCell.cellRow, 2)));

                        double dist_2 = (Math.sqrt(Math.pow(sourceCell.cellCol - previousCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - previousCell.cellRow, 2)));

                        Tuple<Double, Double> intersectionPoint = null;
                        // if radius is between intersectionCell and intersectionCell - 1, we consider these two cells
                        if (radius > dist_2 && radius < dist) {

                            // here compute the intersection point

                            intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                    sourceX, sourceY, radius,
                                    intersectionCell.cellCol, intersectionCell.cellRow,
                                    previousCell.cellCol, previousCell.cellRow);


                        } else if (radius > dist && radius < dist_1) {
                            // if radius is between intersectionCell and intersectionCell + 1 we consider these two cells

                            // here compute the intersection

                            intersectionPoint = computeIntersectionOfCircleAndLineSegment(
                                    sourceX, sourceY, radius,
                                    intersectionCell.cellCol, intersectionCell.cellRow,
                                    nextCell.cellCol, nextCell.cellRow);

                        } else if (radius == dist) {
                            intersectionPoint = new Tuple<Double, Double>((double) intersectionCell.cellCol, (double) intersectionCell.cellRow);
                        }

                        if (intersectionPoint == null) {
                            // System.out.println();
                        }

                        double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionPoint.first, 2) +
                                Math.pow(cell.cellRow - intersectionPoint.second, 2)));

                        double distanceFromSourceToIntersection = (Math.sqrt(Math.pow(sourceCell.cellCol - intersectionPoint.first, 2) +
                                Math.pow(sourceCell.cellRow - intersectionPoint.second, 2)));

                        double angle = Math.acos((Math.pow(radius, 2) + Math.pow(distanceFromSourceToIntersection, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                (2.0 * radius * distanceFromSourceToIntersection));

                        distanceForPathMatrix.distanceMatrix[i][j] = angle;

                    } else {

                        double distanceFromCellToIntersection = (Math.sqrt(Math.pow(cell.cellCol - intersectionCell.cellCol, 2) +
                                Math.pow(cell.cellRow - intersectionCell.cellRow, 2)));

                        double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(distanceFromCellToIntersection, 2)) /
                                (2.0 * radius * radius));

                        distanceForPathMatrix.distanceMatrix[i][j] = angle;

                    }

                }
            }

            Iterator pathCellIter = pathCells.iterator();

            while (pathCellIter.hasNext()) {

                Cell cell = (Cell) pathCellIter.next();

                distanceForPathMatrix.distanceMatrix[cell.cellCol][cell.cellRow] = 0.0;
                //distances[cell.cellCol][cell.cellRow] = 0.0;

            }

            //distancesForPaths.add(distances);
            double[][] transposedMatrix = transposeMatrix(distanceForPathMatrix.distanceMatrix);
            distanceForPathMatrix.distanceMatrix = transposedMatrix;
            distanceForPathMatrix.pathId = path.id;

            distancesForPaths.add(distanceForPathMatrix);
        }
        return distancesForPaths;


    }


    public static Tuple<Double, Double> computeIntersectionOfCircleAndLineSegment(double circlex, double circley, double radius,
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

    public static ArrayList computeVerticalDistances(Cell[][] grid, ArrayList<Path> paths) throws IOException {
        System.out.println("computing vertical distances ");
        logFileWriter.write("computing vertical distances " + "\n");

        Iterator pathIterator = paths.iterator();

        ArrayList distancesForPaths = new ArrayList();

        int pathCounter = 0;

        // for all paths
        while (pathIterator.hasNext()) {

            Path path = (Path) pathIterator.next();

            ArrayList pathCells = path.cells;

            Collections.reverse(pathCells);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    // inefficient method:

                    Iterator pathCellIterator = pathCells.iterator();

                    while (pathCellIterator.hasNext()) {

                        Cell pathCell = (Cell) pathCellIterator.next();

                        if (pathCell.cellCol == cell.cellCol) {

                            double distance = Math.abs(pathCell.cellRow - cell.cellRow);
                            distances[i][j] = distance;
                        }

                    }

                }
            }

            Iterator pathCellIter = pathCells.iterator();

            while (pathCellIter.hasNext()) {

                Cell cell = (Cell) pathCellIter.next();

                distances[cell.cellCol][cell.cellRow] = 0.0;

            }

            double[][] transposedDistancesMatrix = transposeMatrix(distances);
            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distanceMatrix = transposedDistancesMatrix;
            distanceForPathMatrix.pathId = pathCounter;
            pathCounter++;
            distancesForPaths.add(distanceForPathMatrix);
        }
        return distancesForPaths;

    }

    public static ArrayList computeDistancesBruteForce(Cell[][] grid, ArrayList<Path> paths) {

        Iterator<Path> pathIterator = paths.iterator();
        ArrayList distancesForPaths = new ArrayList();

        while (pathIterator.hasNext()) {

            Path path = pathIterator.next();
            ArrayList<Cell> cells = path.cells;

            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.pathId = path.id;
            distanceForPathMatrix.distanceMatrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            for (int col = 0; col < NR_OF_COLUMNS; col++) {
                for (int row = 0; row < NR_OF_ROWS; row++) {

                    Cell cell = grid[col][row];

                    double minDistance = Integer.MAX_VALUE;

                    for (int c = 0; c < cells.size(); c++) {

                        Cell pathCell = cells.get(c);

                        double distance = Math.sqrt(Math.pow(cell.cellCol - pathCell.cellCol, 2) + Math.pow(cell.cellRow - pathCell.cellRow, 2));

                        if (distance < minDistance) {
                            minDistance = distance;
                        }
                    }
                    distanceForPathMatrix.distanceMatrix[col][row] = minDistance;
                }
            }
            double[][] transposed = transposeMatrix(distanceForPathMatrix.distanceMatrix);
            distanceForPathMatrix.distanceMatrix = transposed;
            distancesForPaths.add(distanceForPathMatrix);
        }
        return distancesForPaths;
    }

    public static ArrayList compute_Dijkstra(Cell[][] grid, ArrayList<Path> paths) throws IOException {

        System.out.println("Computing Dijkstra");
        logFileWriter.write("Computing Dijkstra" + "\n");

        Iterator pathIterator = paths.iterator();

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
        ArrayList distancesForPaths = new ArrayList();

        while (pathIterator.hasNext()) {

            Path path = (Path) pathIterator.next();

            // store results for a single path here
            //double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distanceMatrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];
            distanceForPathMatrix.pathId = path.id;

            // Initialize distances with max distance values
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {
                    distanceForPathMatrix.distanceMatrix[j][i] = Integer.MAX_VALUE;
                    grid[j][i].distance = Integer.MAX_VALUE;
                }
            }

            PriorityQueue<Cell> Q = new PriorityQueue<Cell>(NR_OF_COLUMNS * NR_OF_ROWS,
                    new distanceComparator());

            // Loop over each cell on a path:
            Iterator cellIterator = path.cells.iterator();
            while (cellIterator.hasNext()) {

                // A cell on a path:
                Cell cell = (Cell) cellIterator.next();

                //cell.distance = 0.0;
                distanceForPathMatrix.distanceMatrix[cell.cellCol][cell.cellRow] = 0.0;
                cell.distance = 0.0;

                Q.add(cell);

            }

            while (!Q.isEmpty()) {

                //Tuple<Cell, Double> nd = Q.poll();
                Cell cell = Q.poll();

                int x = cell.cellCol;
                int y = cell.cellRow;

                for (int i = 0; i < 16; i++) {

                    int adjX = x + dCol[i];
                    int adjY = y + dRow[i];
                    double weight = weights[i];

                    if (isInsideGrid(adjX, adjY)) {

                        if (distanceForPathMatrix.distanceMatrix[adjX][adjY] > distanceForPathMatrix.distanceMatrix[x][y] + weight) {

                            // If Cell is already been reached once,
                            // remove it from priority queue
                            if (distanceForPathMatrix.distanceMatrix[adjX][adjY] != Integer.MAX_VALUE) {
                                Cell adj = grid[adjX][adjY];//new Cell(rows, cols, dist[rows][cols]);
                                adj.distance = distanceForPathMatrix.distanceMatrix[adjX][adjY];
                                Q.remove(adj);

                            }

                            // Insert cell with updated distance
                            distanceForPathMatrix.distanceMatrix[adjX][adjY] = ((distanceForPathMatrix.distanceMatrix[cell.cellCol][cell.cellRow] + weight));

                            grid[adjX][adjY].distance = distanceForPathMatrix.distanceMatrix[adjX][adjY];
                            Q.add(grid[adjX][adjY]); //new Cell(rows, cols, dist[rows][cols]));

                        }
                    }
                }
            }
            distanceForPathMatrix.distanceMatrix = transposeMatrix(distanceForPathMatrix.distanceMatrix);
            distancesForPaths.add(distanceForPathMatrix);
        }

        return distancesForPaths;
    }

    public static Cell[][] transposeGrid(Cell[][] grid) {
        // TODO: figure out how to assign these without transpose
        int m = grid.length;
        int n = grid[0].length;

        Cell[][] transposedMatrix = new Cell[n][m];

        for (int x = 0; x < n; x++) {
            for (int y = 0; y < m; y++) {
                transposedMatrix[x][y] = grid[y][x];
            }
        }

        return transposedMatrix;
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

    public static ArrayList computeAngularDistance(Cell[][] grid, ArrayList paths) throws IOException {

        System.out.println("computing angular distance");
        logFileWriter.write("computing angular distance" + "\n");

        Iterator pathIterator = paths.iterator();

        ArrayList distancesForPaths = new ArrayList();

        // for all paths
        while (pathIterator.hasNext()) {

            ArrayList path = (ArrayList) pathIterator.next();

            Collections.reverse(path);

            double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            // for each cell in the grid
            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    Cell cell = grid[i][j];

                    double radius = (Math.sqrt(Math.pow(sourceCell.cellCol - cell.cellCol, 2) +
                            Math.pow(sourceCell.cellRow - cell.cellRow, 2)));

                    int indexOfCell = binarySearch_2(path, radius);//binarySearch(path, 0, path.size(), radius);

                    Cell intersectionCell = (Cell) path.get(indexOfCell);

                    // Here we can either use radius as dist or the actual distance. Matter of precision.
                    double dist = (Math.sqrt(Math.pow(intersectionCell.cellCol - sourceCell.cellCol, 2) +
                            Math.pow(intersectionCell.cellRow - sourceCell.cellRow, 2)));

                    double dist_2 = (Math.sqrt(Math.pow(intersectionCell.cellCol - cell.cellCol, 2) +
                            Math.pow(intersectionCell.cellRow - cell.cellRow, 2)));

                    double angle = Math.acos((Math.pow(radius, 2) + Math.pow(dist, 2) - Math.pow(dist_2, 2)) /
                            (2.0 * radius * dist));

                    distances[i][j] = angle;
                }
            }

            Iterator pathCellIter = path.iterator();

            while (pathCellIter.hasNext()) {

                Cell cell = (Cell) pathCellIter.next();

                distances[cell.cellCol][cell.cellRow] = 0.0;

            }

            distancesForPaths.add(transposeMatrix(distances));
        }
        return distancesForPaths;

    }

    public static int binarySearch_2(ArrayList<Cell> path, double distanceTarget) {

        int start = 0;
        int end = path.size() - 1;

        int ans = -1;

        while (start <= end) {
            int mid = (start + end) / 2;
            Cell midCell = path.get(mid);

            double distToMid = (Math.sqrt(Math.pow(sourceCell.cellCol - midCell.cellCol, 2) +
                    Math.pow(sourceCell.cellRow - midCell.cellRow, 2)));

            if (distToMid == distanceTarget) {
                return mid;
            }

            if (distToMid <= distanceTarget) {
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

            Cell midCell = path.get(mid);

            double distToMid = (Math.sqrt(Math.pow(sourceCell.cellCol - midCell.cellCol, 2) +
                    Math.pow(sourceCell.cellRow - midCell.cellRow, 2)));

            if (distToMid == dist) {
                return mid;
            } else if (distToMid > dist) {
                binarySearch(path, l, mid - 1, dist);
            } else {
                binarySearch(path, mid + 1, r, dist);
            }
        }

        return -1;
    }

    public static ArrayList computeAngularDistance_2(Cell[][] grid, ArrayList paths) {

        System.out.println("computing angular distance");

        Iterator pathIterator = paths.iterator();

        ArrayList distancesForPaths = new ArrayList();

        Cell sourceCell = grid[sourceX][sourceY];

        // for all paths
        while (pathIterator.hasNext()) {

            ArrayList path = (ArrayList) pathIterator.next();

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

                    double radius = (Math.sqrt(Math.pow(sourceCell.cellCol - cell.cellCol, 2) +
                            Math.pow(sourceCell.cellRow - cell.cellRow, 2)));

                    Iterator cellIterator = path.iterator();

                    while (cellIterator.hasNext()) {

                        Cell pathCell = (Cell) cellIterator.next();

                        double dist = (Math.sqrt(Math.pow(sourceCell.cellCol - pathCell.cellCol, 2) +
                                Math.pow(sourceCell.cellRow - pathCell.cellRow, 2)));
                        if (dist <= radius) {
                            // consider next cell on path
                            continue;
                        } else {
                            // this is the cell
                            double dist_2 = (Math.sqrt(Math.pow(cell.cellCol - pathCell.cellCol, 2) +
                                    Math.pow(cell.cellRow - pathCell.cellRow, 2)));

                            double angle = Math.acos((Math.pow(dist, 2) + Math.pow(radius, 2) - Math.pow(dist_2, 2)) /
                                    (2.0 * dist * radius));

                            distances[i][j] = angle;

                        }
                    }
                }
            }

            distancesForPaths.add(transposeMatrix(distances));
        }
        return distancesForPaths;
    }

    public static ArrayList computeBfs(Cell[][] grid, ArrayList<Path> paths) throws IOException {

        System.out.println("computing bfs");
        logFileWriter.write("computing bfs" + "\n");

        Iterator pathIterator = paths.iterator();

        // Direction vectors
        int dRow[] = {-1, 0, 1, 0};
        int dCol[] = {0, 1, 0, -1};

        // dw = weight of the edge

        ArrayList distancesForPaths = new ArrayList();

        while (pathIterator.hasNext()) {


            boolean[][] visited = new boolean[NR_OF_COLUMNS][NR_OF_ROWS];

            //double[][] distances = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            DistanceForPathMatrix distanceForPathMatrix = new DistanceForPathMatrix();
            distanceForPathMatrix.distanceMatrix = new double[NR_OF_COLUMNS][NR_OF_ROWS];

            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    distanceForPathMatrix.distanceMatrix[i][j] = 0.0f;

                }
            }

            Path path = (Path) pathIterator.next();

            distanceForPathMatrix.pathId = path.id;

            Queue<Cell> queue = new LinkedList<>();

            Iterator cellIterator = path.cells.iterator();

            while (cellIterator.hasNext()) {

                Cell cell = (Cell) cellIterator.next();

                // add all cells of a path to queue
                queue.add(cell);

                visited[cell.cellRow][cell.cellCol] = true;

            }

            while (!queue.isEmpty()) {

                Cell cell = queue.peek();

                int x = cell.cellCol;
                int y = cell.cellRow;

                queue.remove();

                for (int i = 0; i < 4; i++) {

                    int adjX = x + dCol[i];
                    int adjY = y + dRow[i];

                    if (isValid(visited, adjY, adjX)) {

                        queue.add(grid[adjX][adjY]);
                        visited[adjY][adjX] = true;

                        distanceForPathMatrix.distanceMatrix[adjY][adjX] = distanceForPathMatrix.distanceMatrix[y][x] + 1;
                    }
                }
            }

            distancesForPaths.add(distanceForPathMatrix);
        }

        return distancesForPaths;
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

    public static void computeShortestPathsNaive(Cell[][] grid, ArrayList paths) {

        //https://programmer.ink/think/graph-theory-search-how-to-use-multi-source-bfs-to-reduce-time-complexity.html
        //https://www.geeksforgeeks.org/multi-source-shortest-path-in-unweighted-graph/

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                Cell currentCell = grid[i][j];

                Iterator pathIt = paths.iterator();

                while (pathIt.hasNext()) {

                    ArrayList path = (ArrayList) pathIt.next();

                    Iterator pathIter = path.iterator();

                    ArrayList distances = new ArrayList();

                    while (pathIter.hasNext()) {

                        Cell pathCell = (Cell) pathIter.next();

                        double distance = (double) (Math.sqrt(Math.pow(pathCell.cellCol - currentCell.cellCol, 2) + Math.pow(pathCell.cellRow - currentCell.cellRow, 2)));

                        distances.add(distance);

                    }

                    int minDistanceIndex = distances.indexOf(Collections.min(distances));

                    Cell minDistanceCell = (Cell) path.get(minDistanceIndex);

                }
            }
        }
    }

    public static void drawDistances(Cell[][] grid, ArrayList paths, ArrayList distancesForPaths,
                                     boolean showIntermediateResults, double width, double scale, int imageIndex,
                                     String iterationLocation) throws
            IOException {

//        jframe = new JFrame("panel");
//        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);
//
//        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
//                BufferedImage.TYPE_INT_ARGB);
//
//        double[][] dist = (double[][]) distancesForPaths.get(1);
//
//        double maxDist = dist[0][0];
//        double minDist = 0;
//
//        for (int i = 0; i < NR_OF_COLUMNS; i++) {
//
//            for (int j = 0; j < NR_OF_ROWS; j++) {
//
//                if (dist[i][j] > maxDist) {
//                    maxDist = dist[i][j];
//                }
//            }
//        }
//
//        for (int i = 0; i < NR_OF_COLUMNS; i++) {
//            for (int j = 0; j < NR_OF_ROWS; j++) {
//
//                int c = (int) (dist[i][j] * 255.0 / maxDist);
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
//            Iterator cellIter = path.iterator();
//
//            while (cellIter.hasNext()) {
//
//                Cell cell = (Cell) cellIter.next();
//
//                image.setRGB((int) cell.cellX, (int) cell.cellY, new Color(255, 255, 255).getRGB());
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

        double[][] dist = (double[][]) distancesForPaths.get(0);

        double maxDist = dist[0][0];
        double minDist = 0;

        double[][] minDistances = new double[NR_OF_ROWS][NR_OF_COLUMNS];

        for (int i = 0; i < NR_OF_ROWS; i++) {
            for (int j = 0; j < NR_OF_COLUMNS; j++) {

                Iterator distIter = distancesForPaths.iterator();

                double minDistance = (double) dist[i][j];

                while (distIter.hasNext()) {

                    double[][] distancesForPath = (double[][]) distIter.next();

                    if (distancesForPath[i][j] < minDistance) {

                        minDistance = distancesForPath[i][j];

                    }

                }
                minDistances[i][j] = minDistance;

            }
        }

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                if (minDistances[i][j] > maxDist) {
                    maxDist = minDistances[i][j];
                }
            }
        }

        System.out.println("min dist: " + minDist + " max dist : " + maxDist);
        logFileWriter.write("min dist: " + minDist + " max dist : " + maxDist + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                //int c = (int) (dist[i][j] * 255.0 / maxDist);

                float value = (float) ((minDistances[j][i] - minDist) / (maxDist - minDist));

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

                Iterator cellIter = path.cells.iterator();

                while (cellIter.hasNext()) {

                    Cell cell = (Cell) cellIter.next();

                    image.setRGB((int) cell.cellCol, (int) cell.cellRow, new Color(0, 0, 0).getRGB());

                }
            }
        }

        // draw string in image:

        if (DRAW_TEXT_DESCRIPTION) {
            Font f = new Font(Font.MONOSPACED, Font.PLAIN, 20);
            //String s = "width: " + width + " scale: " + scale + " i: " + imageIndex;
            String s = "width: " + width + " i: " + imageIndex;
            Graphics g = image.getGraphics();
            g.setColor(Color.BLUE);
            g.setFont(f);
            FontMetrics fm = g.getFontMetrics();
            int x = image.getWidth() - fm.stringWidth(s) - 5;
            int y = fm.getHeight();
            g.drawString(s, x, y);
            g.dispose();
        }

        if (showIntermediateResults) {

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
        //dir = new File(currentWorkingPath.concat("\\" +  storageLocationName + "\\" + iterationLocation + "\\"));

        File file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/imageDistances_" + imageIndex + ".png"));
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

        Color resultColor = new Color(red, green, blue);

        return resultColor;
    }

    public static void drawObstacles(double[][] matrix, ArrayList paths, int imageIndex, String iterationLocation) throws IOException {
        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_COLUMNS, NR_OF_ROWS);

        BufferedImage image = new BufferedImage(NR_OF_COLUMNS, NR_OF_ROWS,
                BufferedImage.TYPE_INT_ARGB);

        double maxHeight = matrix[0][0];
        double minHieght = matrix[0][0];

        maxHeight = MAX_HEIGHT;
        minHieght = MIN_HEIGHT;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                if ( HEIGHT_UPDATE[i][j] + matrix[i][j] < minHieght) {
                    minHieght = HEIGHT_UPDATE[i][j] +  matrix[i][j];
                }

                if ( HEIGHT_UPDATE[i][j] + matrix[i][j] > maxHeight) {
                    maxHeight = HEIGHT_UPDATE[i][j] +  matrix[i][j];
                }

            }
        }

        System.out.println("min height update : " + minHieght + " max height update : " + maxHeight);
        logFileWriter.write("min height update : " + minHieght + " max height update : " + maxHeight + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                double heightUpdate = HEIGHT_UPDATE[i][j];

                float value = (float) ((heightUpdate + matrix[i][j] - minHieght) / (maxHeight - minHieght));

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
                image.setRGB(i, NR_OF_ROWS - 1 - j, color.getRGB());

                //  } catch (Exception e) {
                //      System.out.println();
                //  }

            }
        }


        Iterator iter = paths.iterator();

        while (iter.hasNext()) {

            Path path = (Path) iter.next();

            Iterator cellIter = path.cells.iterator();

            while (cellIter.hasNext()) {

                Cell cell = (Cell) cellIter.next();

                image.setRGB((int) cell.cellCol, NR_OF_ROWS - 1 - (int) cell.cellRow, new Color(0, 0, 0).getRGB());

            }
        }

        File file;
        file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/obstacles_" + imageIndex + ".png"));

        file.mkdirs();
        ImageIO.write(image, "png", file);
    }

    public static void drawMatrix(double[][] matrix, ArrayList paths, int imageIndex, String iterationLocation, boolean relativeToTotal) throws IOException {
        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        double maxHeight = matrix[0][0];
        double minHieght = matrix[0][0];

        if (relativeToTotal) {
            maxHeight = MAX_HEIGHT;
            minHieght = MIN_HEIGHT;

            for (int i = 0; i < NR_OF_COLUMNS; i++) {
                for (int j = 0; j < NR_OF_ROWS; j++) {

                    if (matrix[i][j] < minHieght) {
                        minHieght = matrix[i][j];
                    }

                    if (matrix[i][j] > maxHeight) {
                        maxHeight = matrix[i][j];
                    }

                }
            }

        } else {
            for (int i = 0; i < NR_OF_COLUMNS; i++) {

                for (int j = 0; j < NR_OF_ROWS; j++) {

                    if (matrix[i][j] > maxHeight) {
                        maxHeight = matrix[i][j];
                    }
                    if (matrix[i][j] < minHieght) {
                        minHieght = matrix[i][j];
                    }
                }
            }
            double temp = minHieght;
            minHieght = maxHeight;
            maxHeight = temp;
        }

        System.out.println("min height update : " + minHieght + " max height update : " + maxHeight);
        logFileWriter.write("min height update : " + minHieght + " max height update : " + maxHeight + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                float value = (float) ((matrix[i][j] - minHieght) / (maxHeight - minHieght));

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

        if (DRAW_PATHS) {

            Iterator iter = paths.iterator();

            while (iter.hasNext()) {

                Path path = (Path) iter.next();

                Iterator cellIter = path.cells.iterator();

                while (cellIter.hasNext()) {

                    Cell cell = (Cell) cellIter.next();

                    image.setRGB((int) cell.cellCol, (int) cell.cellRow, new Color(0, 0, 0).getRGB());

                }
            }
        }
        File file;
        if (relativeToTotal) {
            file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/updateGlobalHeight_" + imageIndex + ".png"));
        } else {
            file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/updateLocalHeight_" + imageIndex + ".png"));
        }
        file.mkdirs();
        ImageIO.write(image, "png", file);
    }

    public static void drawFlowAccumulation(Cell[][] grid, ArrayList paths, int imageIndex, boolean showIntermediateResults, String iterationLocation, double width, double scale)
            throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        double maxHeight = grid[0][0].flowAccumulation;
        double minHieght = grid[0][0].flowAccumulation;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                if (grid[i][j].flowAccumulation > maxHeight) {
                    maxHeight = grid[i][j].flowAccumulation;
                }
                if (grid[i][j].flowAccumulation < minHieght) {
                    minHieght = grid[i][j].flowAccumulation;
                }
            }
        }

        System.out.println("min flowAccumulation : " + minHieght + " max flowAccumulation : " + maxHeight);
        logFileWriter.write("min flowAccumulation : " + minHieght + " max flowAccumulation : " + maxHeight + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                float value = (float) ((grid[i][j].flowAccumulation - minHieght) / (maxHeight - minHieght));

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

                image.setRGB(i, j, color.getRGB());

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
//                Iterator cellIter = path.cells.iterator();
//
//                while (cellIter.hasNext()) {
//
//                    Cell cell = (Cell) cellIter.next();
//
//                    image.setRGB((int) cell.cellCol, (int) cell.cellRow, new Color(0, 0, 0).getRGB());
//
//                }
//            }
//        }


        // draw string in image:

        if (DRAW_TEXT_DESCRIPTION) {
            Font f = new Font(Font.MONOSPACED, Font.PLAIN, 20);
            String s = "width: " + width + " scale: " + scale + " i: " + imageIndex;
            Graphics g = image.getGraphics();
            g.setColor(Color.BLUE);
            g.setFont(f);
            FontMetrics fm = g.getFontMetrics();
            int x = image.getWidth() - fm.stringWidth(s) - 5;
            int y = fm.getHeight();
            g.drawString(s, x, y);
            g.dispose();
        }

        if (showIntermediateResults) {

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
        //dir = new File(currentWorkingPath.concat("\\" +  storageLocationName + "\\" + iterationLocation + "\\"));

        File file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/flowAccumulation_" + imageIndex + ".png"));
        file.mkdirs();
        ImageIO.write(image, "png", file);

    }

    public static void drawFlow(Cell[][] grid, ArrayList paths, int imageIndex, boolean showIntermediateResults, String iterationLocation, double width, double scale)
            throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                int flow = grid[i][j].flowDirection;
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

        if (showIntermediateResults) {

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
        //dir = new File(currentWorkingPath.concat("\\" +  storageLocationName + "\\" + iterationLocation + "\\"));

        File file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/paths_" + imageIndex + ".png"));
        file.mkdirs();
        ImageIO.write(image, "png", file);

    }

    public static void drawPaths(Cell[][] grid, ArrayList paths, int imageIndex, boolean showIntermediateResults, String iterationLocation, double width, double scale)
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

                Iterator cellIter = path.cells.iterator();

                while (cellIter.hasNext()) {

                    Cell cell = (Cell) cellIter.next();

                    image.setRGB((int) cell.cellCol, (int) cell.cellRow, new Color(255, 255, 255).getRGB());

                }
            }
        }

        if (showIntermediateResults) {

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
        //dir = new File(currentWorkingPath.concat("\\" +  storageLocationName + "\\" + iterationLocation + "\\"));

        File file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/paths_" + imageIndex + ".png"));
        file.mkdirs();
        ImageIO.write(image, "png", file);

    }

    public static void draw(Cell[][] grid, ArrayList paths, int imageIndex, boolean showIntermediateResults, String iterationLocation, double width, double scale, ArrayList<Point> pointsList)
            throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(NR_OF_ROWS, NR_OF_COLUMNS);

        BufferedImage image = new BufferedImage(NR_OF_ROWS, NR_OF_COLUMNS,
                BufferedImage.TYPE_INT_ARGB);

        double maxHeight = grid[0][0].height;
        double minHieght = grid[0][0].height;

        for (int i = 0; i < NR_OF_COLUMNS; i++) {

            for (int j = 0; j < NR_OF_ROWS; j++) {

                if (grid[i][j].height > maxHeight) {
                    maxHeight = grid[i][j].height;
                }
                if (grid[i][j].height < minHieght) {
                    minHieght = grid[i][j].height;
                }
            }
        }

        System.out.println("min height : " + minHieght + " max height : " + maxHeight);
        logFileWriter.write("min height : " + minHieght + " max height : " + maxHeight + "\n");

        for (int i = 0; i < NR_OF_COLUMNS; i++) {
            for (int j = 0; j < NR_OF_ROWS; j++) {

                float value = (float) ((grid[i][j].height - minHieght) / (maxHeight - minHieght));

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

                image.setRGB(i, j, color.getRGB());

            }
        }

        if (DRAW_PATHS) {

            Iterator iter = paths.iterator();

            while (iter.hasNext()) {

                Path path = (Path) iter.next();

                Iterator cellIter = path.cells.iterator();

                while (cellIter.hasNext()) {

                    Cell cell = (Cell) cellIter.next();

                    image.setRGB((int) cell.cellCol, (int) cell.cellRow, new Color(0, 0, 0).getRGB());

                }
            }
            for (int i = 0; i < pointsList.size(); i++) {

                if (pointsList.get(i).name.equals("S")) {
                    continue;
                }

                image.setRGB((int) pointsList.get(i).gridCol, (int) pointsList.get(i).gridRow, new Color(255, 255, 255).getRGB());

                if (pointsList.get(i).gridCol - 1 >= 0) {
                    image.setRGB(pointsList.get(i).gridCol - 1, (int) pointsList.get(i).gridRow, new Color(255, 255, 255).getRGB());
                }
                if (pointsList.get(i).gridCol + 1 <= NR_OF_COLUMNS - 1) {
                    image.setRGB(pointsList.get(i).gridCol + 1, (int) pointsList.get(i).gridRow, new Color(255, 255, 255).getRGB());
                }
                if ((int) pointsList.get(i).gridRow - 1 >= 0) {
                    image.setRGB(pointsList.get(i).gridCol, (int) pointsList.get(i).gridRow - 1, new Color(255, 255, 255).getRGB());
                }
                if ((int) pointsList.get(i).gridRow + 1 <= NR_OF_ROWS - 1) {
                    image.setRGB(pointsList.get(i).gridCol, (int) pointsList.get(i).gridRow + 1, new Color(255, 255, 255).getRGB());
                }

            }
        }


        // draw string in image:
        if (DRAW_TEXT_DESCRIPTION) {
            Font f = new Font(Font.MONOSPACED, Font.PLAIN, 20);
            //String s = "width: " + width + " scale: " + scale + " i: " + imageIndex;
            String s = "width: " + width + " i: " + imageIndex;
            Graphics g = image.getGraphics();
            g.setColor(Color.BLUE);
            g.setFont(f);
            FontMetrics fm = g.getFontMetrics();
            int x = image.getWidth() - fm.stringWidth(s) - 5;
            int y = fm.getHeight();
            g.drawString(s, x, y);
            g.dispose();
        }

        if (showIntermediateResults) {

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
        //dir = new File(currentWorkingPath.concat("\\" +  storageLocationName + "\\" + iterationLocation + "\\"));

        int w = (int) width;
        int s = (int) scale;

        File file = new File(currentWorkingPath.concat("/" + storageLocationName + "/" + iterationLocation + "/height_w_" + w + "_s_" + s + "_" + imageIndex + ".png"));
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

    public static ArrayList computePathsToFrameEdge(ArrayList pointsList, Cell[][] grid) throws IOException {
        System.out.println("computing paths from fram to edge");
        logFileWriter.write("computing paths from fram to edge" + "\n");

        Iterator it = pointsList.iterator();

        // ArrayList<ArrayList<Cell>> paths = new ArrayList();
        ArrayList<Path> paths = new ArrayList();

        int pathIdCounter = 0;
        while (it.hasNext()) {

            Point point = (Point) it.next();

            if (point.name.equals(TARGET_NAME)) {
                continue;
            }

            Cell gridCell = grid[point.gridCol][point.gridRow];

            Cell currentCell = gridCell;

            ArrayList<Cell> path = new ArrayList();

            path.add(currentCell);

            int counter = 0;

            while (!(currentCell.title.equals("rightEdge"))) {

                if (counter > 4 * (NR_OF_COLUMNS + NR_OF_ROWS)) {
                    System.out.println("something went wrong");
                    logFileWriter.write("something went wrong" + "\n");
                    return null;
                }

                int col = (int) currentCell.cellCol;
                int row = (int) currentCell.cellRow;

                double flow = currentCell.flowDirection;

                if (flow == 1) {
                    currentCell = grid[col + 1][row];
                } else if (flow == 2) {

                    if (row + 1 == NR_OF_ROWS && CIRCULAR_MODE) {
                        row = -1;
                    }

                    currentCell = grid[col + 1][row + 1];
                } else if (flow == 4) {

                    if (row + 1 == NR_OF_ROWS && CIRCULAR_MODE) {
                        row = -1;
                    }

                    currentCell = grid[col][row + 1];
                } else if (flow == 8) {

                    if (row + 1 == NR_OF_ROWS && CIRCULAR_MODE) {
                        row = -1;
                    }

                    currentCell = grid[col - 1][row + 1];
                } else if (flow == 16) {
                    currentCell = grid[col - 1][row];
                } else if (flow == 32) {

                    if (row - 1 == -1 && CIRCULAR_MODE) {
                        row = 500;
                    }

                    currentCell = grid[col - 1][row - 1];
                } else if (flow == 64) {

                    if (row - 1 == -1 && CIRCULAR_MODE) {
                        row = 500;
                    }

                    currentCell = grid[col][row - 1];
                } else if (flow == 128) {

                    if (row - 1 == -1 && CIRCULAR_MODE) {
                        row = 500;
                    }

                    currentCell = grid[col + 1][row - 1];
                }

                path.add(currentCell);

                counter++;

            }

            Path newPath = new Path();
            newPath.cells = path;
            newPath.id = pathIdCounter;
            paths.add(newPath);
            pathIdCounter++;
        }
        return paths;
    }

    public static ArrayList<Path> computePaths(ArrayList pointsList, Cell[][] grid) throws IOException {

        System.out.println("computing paths");
        logFileWriter.write("computing paths" + "\n");

        Iterator it = pointsList.iterator();

        ArrayList<Path> paths = new ArrayList();

        int idCounter = 0;
        while (it.hasNext()) {

            Point point = (Point) it.next();

            if (point.name.equals(TARGET_NAME) || point.name.equals("S")) {
                continue;
            }

            Cell gridCell = grid[point.gridCol][point.gridRow];

            Cell currentCell = gridCell;

            Path path = new Path();
            path.cells = new ArrayList();
            path.id = idCounter;
            //ArrayList<Cell> path = new ArrayList();

            path.cells.add(currentCell);

            int counter = 0;

            while (!(currentCell.title.equals(TARGET_NAME))) {

                if (counter > 4 * (NR_OF_COLUMNS + NR_OF_ROWS)) {
                    System.out.println("something went wrong");
                    logFileWriter.write("something went wrong" + "\n");

                    return null;
                }

                int col = (int) currentCell.cellCol;
                int row = (int) currentCell.cellRow;

                double flow = currentCell.flowDirection;

                if (flow == 1) {
                    currentCell = grid[col + 1][row];
                } else if (flow == 2) {

                    if (row + 1 == NR_OF_ROWS && CIRCULAR_MODE) {
                        row = -1;
                    }

                    currentCell = grid[col + 1][row + 1];
                } else if (flow == 4) {

                    if (row + 1 == NR_OF_ROWS && CIRCULAR_MODE) {
                        row = -1;
                    }

                    currentCell = grid[col][row + 1];
                } else if (flow == 8) {

                    if (row + 1 == NR_OF_ROWS && CIRCULAR_MODE) {
                        row = -1;
                    }

                    currentCell = grid[col - 1][row + 1];
                } else if (flow == 16) {
                    currentCell = grid[col - 1][row];
                } else if (flow == 32) {

                    if (row - 1 == -1 && CIRCULAR_MODE) {
                        row = 500;
                    }

                    currentCell = grid[col - 1][row - 1];
                } else if (flow == 64) {

                    if (row - 1 == -1 && CIRCULAR_MODE) {
                        row = 500;
                    }

                    currentCell = grid[col][row - 1];
                } else if (flow == 128) {

                    if (row - 1 == -1 && CIRCULAR_MODE) {
                        row = 500;
                    }

                    currentCell = grid[col + 1][row - 1];
                }

                path.cells.add(currentCell);

                counter++;

            }
            idCounter++;

            paths.add(path);

        }
        return paths;
    }

    public static void computeFlow(Cell[][] grid, int iteration) throws IOException {
        System.out.println("computing flow");
        logFileWriter.write("computing flow" + "\n");

        for (int col = 0; col < NR_OF_COLUMNS; col++) {

            for (int row = 0; row < NR_OF_ROWS; row++) {

//                int col = (int) grid[col][row].cellCol;
//                int row = (int) grid[col][row].cellRow;

                ArrayList<Cell> neighbors = new ArrayList();

                if (col == sourceX && row == sourceY) {
                    continue;
                }

                if (grid[col][row].title.equals("rightEdge") && HORIZONTAL_FLOW_MODE) {
                    continue;
                }

                Cell left = null;
                Cell right = null;
                Cell top = null;
                Cell bottom = null;
                Cell topLeft = null;
                Cell topRight = null;
                Cell bottomLeft = null;
                Cell bottomRight = null;

                if (CIRCULAR_MODE) {

                    if (row == 0) {
                        top = grid[col][NR_OF_ROWS - 1];
                        neighbors.add(top);
                    }

                    if (row == NR_OF_ROWS - 1) {
                        bottom = grid[col][0];
                        neighbors.add(bottom);
                    }

                    if (row == 0 && col != NR_OF_COLUMNS - 1) {
                        topRight = grid[col + 1][NR_OF_ROWS - 1];
                        neighbors.add(topRight);
                    }

                    if (row == 0 && col != 0) {
                        topLeft = grid[col - 1][NR_OF_ROWS - 1];
                        neighbors.add(topLeft);
                    }

                    if (row == NR_OF_ROWS - 1 && col != NR_OF_COLUMNS - 1) {
                        bottomRight = grid[col + 1][0];
                        neighbors.add(bottomRight);
                    }

                    if (row == NR_OF_ROWS - 1 && col != 0) {
                        bottomLeft = grid[col - 1][0];
                        neighbors.add(bottomLeft);
                    }

                    if (col - 1 >= 0) {
                        left = grid[col - 1][row];
                        neighbors.add(left);
                    }

                    if (row + 1 < NR_OF_ROWS) {
                        bottom = grid[col][row + 1];
                        neighbors.add(bottom);
                    }

                    if (col + 1 < NR_OF_COLUMNS) {
                        right = grid[col + 1][row];
                        neighbors.add(right);
                    }

                    if (row - 1 >= 0) {
                        top = grid[col][row - 1];
                        neighbors.add(top);
                    }

                    if (col - 1 >= 0 && row - 1 >= 0) {
                        topLeft = grid[col - 1][row - 1];
                        neighbors.add(topLeft);
                    }

                    if (row - 1 >= 0 && col + 1 < NR_OF_COLUMNS) {
                        topRight = grid[col + 1][row - 1];
                        neighbors.add(topRight);
                    }

                    if (col - 1 >= 0 && row + 1 < NR_OF_ROWS) {
                        bottomLeft = grid[col - 1][row + 1];
                        neighbors.add(bottomLeft);
                    }

                    if (col + 1 < NR_OF_COLUMNS && row + 1 < NR_OF_ROWS) {
                        bottomRight = grid[col + 1][row + 1];
                        neighbors.add(bottomRight);
                    }

                } else {

                    if (col - 1 >= 0) {
                        left = grid[col - 1][row];
                        neighbors.add(left);
                    }

                    if (row + 1 < NR_OF_ROWS) {
                        bottom = grid[col][row + 1];
                        neighbors.add(bottom);
                    }

                    if (col + 1 < NR_OF_COLUMNS) {
                        right = grid[col + 1][row];
                        neighbors.add(right);
                    }

                    if (row - 1 >= 0) {
                        top = grid[col][row - 1];
                        neighbors.add(top);
                    }

                    if (col - 1 >= 0 && row - 1 >= 0) {
                        topLeft = grid[col - 1][row - 1];
                        neighbors.add(topLeft);
                    }

                    if (row - 1 >= 0 && col + 1 < NR_OF_COLUMNS) {
                        topRight = grid[col + 1][row - 1];
                        neighbors.add(topRight);
                    }

                    if (col - 1 >= 0 && row + 1 < NR_OF_ROWS) {
                        bottomLeft = grid[col - 1][row + 1];
                        neighbors.add(bottomLeft);
                    }

                    if (col + 1 < NR_OF_COLUMNS && row + 1 < NR_OF_ROWS) {
                        bottomRight = grid[col + 1][row + 1];
                        neighbors.add(bottomRight);
                    }
                }

                Iterator it = neighbors.iterator();

                ArrayList dropForNeighbors = new ArrayList();

                while (it.hasNext()) {

                    Cell neighbor = (Cell) it.next();

                    double changeInHeight = grid[col][row].height - neighbor.height;

                    double distance = 0.0;

                    if (neighbor == left || neighbor == right || neighbor == top || neighbor == bottom) {

                        distance = 1.0;

                    } else if (neighbor == topLeft || neighbor == topRight || neighbor == bottomLeft || neighbor == bottomRight) {

                        double distFromSourceToCell = Math.sqrt(Math.pow(grid[row][col].cellCol - sourceX, 2) + Math.pow(grid[row][col].cellRow - sourceY, 2));

                        if (true) {//distFromSourceToCell < 100) {

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

                    double drop = (changeInHeight / distance);
                    dropForNeighbors.add(drop);

                }

                int maxDropIndex = dropForNeighbors.indexOf(Collections.max(dropForNeighbors));

                Cell maxDropNeighbor = (Cell) neighbors.get(maxDropIndex);

                if (maxDropNeighbor == left) {
                    grid[col][row].flowDirection = 16;
                } else if (maxDropNeighbor == right) {
                    grid[col][row].flowDirection = 1;
                } else if (maxDropNeighbor == bottom) {
                    grid[col][row].flowDirection = 4;
                } else if (maxDropNeighbor == top) {
                    grid[col][row].flowDirection = 64;
                } else if (maxDropNeighbor == topLeft) {
                    grid[col][row].flowDirection = 32;
                } else if (maxDropNeighbor == topRight) {
                    grid[col][row].flowDirection = 128;
                } else if (maxDropNeighbor == bottomLeft) {
                    grid[col][row].flowDirection = 8;
                } else if (maxDropNeighbor == bottomRight) {
                    grid[col][row].flowDirection = 2;
                }

            }
        }
    }


    public static void initializePointsInGrid(Cell[][] grid, ArrayList<Point> pointsList) {

        Iterator it = pointsList.iterator();

        // setting the names of source/targets and setting height of source
        while (it.hasNext()) {

            Point point = (Point) it.next();

            grid[point.gridCol][point.gridRow].title = point.name;

            if (point.name.equals(TARGET_NAME)) {
                grid[point.gridCol][point.gridRow].height = -10000;

                sourceX = point.gridCol;
                sourceY = point.gridRow;

            }
        }

        sourceCell = grid[sourceX][sourceY];
    }


    public static void computeCellForPoint(Bounds bounds, ArrayList points) {

        Iterator it = points.iterator();

        double extentX = bounds.maxX - bounds.minX;
        double extentY = bounds.maxY - bounds.minY;

        double stepX = extentX / NR_OF_COLUMNS;
        double stepY = extentY / NR_OF_ROWS;

        while (it.hasNext()) {

            Point point = (Point) it.next();

            int i = 1;
            int j = 1;

            while (bounds.minX + i * stepX < point.col) {
                i = i + 1;
            }
            int columnIndex = i - 1;

            while (bounds.minY + j * stepY < point.row) {
                j = j + 1;
            }
            int rowIndex = j - 1;

            point.gridCol = columnIndex;
            point.gridRow = rowIndex;

        }

    }

    public static ArrayList processInput(ArrayList items) {

        Iterator it = items.iterator();

        ArrayList pointList = new ArrayList();

        while (it.hasNext()) {

            String inputLine = it.next().toString();
            List<String> components = Arrays.asList(inputLine.split(";"));

            Point point = new Point();
            point.name = components.get(0);
            point.col = Float.parseFloat(components.get(1));
            point.row = -Float.parseFloat(components.get(2));

            pointList.add(point);

        }

        return pointList;
    }

    public static Bounds obtainBounds(ArrayList<Point> inputPoints) {

        Iterator it = inputPoints.iterator();

        double minX = inputPoints.get(0).col;
        double minY = inputPoints.get(0).row;
        double maxX = inputPoints.get(0).col;
        double maxY = inputPoints.get(0).row;

        while (it.hasNext()) {

            Point nextPoint = (Point) it.next();

            if (nextPoint.col < minX) {
                minX = nextPoint.col;
            }
            if (nextPoint.col > maxX) {
                maxX = nextPoint.col;
            }
            if (nextPoint.row < minY) {
                minY = nextPoint.row;
            }
            if (nextPoint.row > maxY) {
                maxY = nextPoint.row;
            }

        }

        Bounds result = new Bounds();
        result.maxX = maxX;
        result.minX = minX;
        result.minY = minY;
        result.maxY = maxY;

        return result;
    }


    public static ArrayList readInput() throws FileNotFoundException {

        Scanner sc = new Scanner(new File(INPUT_FILE_NAME));
        sc.useDelimiter("\n");

        ArrayList items = new ArrayList();

        while (sc.hasNext()) {

            items.add(sc.next());

        }
        sc.close();

        return items;
    }


    public static Cell[][] initializeGrid(int nrOfRows, int nrOfColumns) {

        Cell[][] grid = new Cell[nrOfRows][nrOfColumns];

        for (int i = 0; i < nrOfRows; i++) {
            for (int j = 0; j < nrOfColumns; j++) {

                Cell cell = new Cell();

                cell.cellCol = i;
                cell.cellRow = j;
                cell.flowDirection = 0;
                cell.height = 0;

                cell.title = "";
                grid[i][j] = cell;
                cell.isObstacle = false;

                // i = columns
                // j = rows
                if (i == nrOfRows - 1 || j == nrOfColumns - 1 || i == 0 || j == 0) {
                    cell.title = "edge";
                }

                if (i == nrOfRows - 1) {
                    cell.title = "rightEdge";
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
