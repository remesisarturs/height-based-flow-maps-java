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


    public static int nr_of_rows = 500;
    public static int nr_of_columns = 500;

    public static int source_x;
    public static int source_y;

    public static String target_name = "A";

    static JFrame jframe;


    private static double zoomFactor = 1;
    private static double prevZoomFactor = 1;
    private static boolean zoomer;


    public static DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm_ss");
    public static LocalDateTime now = LocalDateTime.now();
    public static String storage_location_name;
    public static String currentWorkingPath;


    public static double height_function_width;
    public static double height_function_scale;

    public static int NR_OF_ITERATIONS;

    public static String BASE_HEIGHT_TYPE;

    public static double base_scale;

    public static void main(String[] args) throws IOException {

        BASE_HEIGHT_TYPE = "SQRT";
        BASE_HEIGHT_TYPE = "OTHER";

        storage_location_name = dtf.format(now);

        currentWorkingPath = System.getProperty("user.dir").concat("\\experiments\\");

        storage_location_name = storage_location_name.concat("_" + BASE_HEIGHT_TYPE);

        File dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\"));

        dir.mkdir();

        Cell[][] grid = initialize_grid(nr_of_rows, nr_of_columns);

        ArrayList<String> items = read_input();

        ArrayList<Point> points_list = process_input(items);

        Bounds bounds = obtain_bounds(points_list);

        compute_cell_for_point(bounds, points_list);

        initialize_points_in_grid(grid, points_list);

        base_scale = 0.5;

        if (BASE_HEIGHT_TYPE.equals("SQRT")) {
            initialize_grid_height_sqrt(grid, base_scale);
        } else {
            initialize_grid_height(grid, base_scale);
        }

        compute_flow(grid);

        ArrayList<ArrayList<Cell>> paths = compute_paths(points_list, grid);

        draw(grid, paths, 0, false);

        ArrayList<double[][]> distances_for_paths = compute_bfs(grid, paths);

        NR_OF_ITERATIONS = 10;

        Tuple<Cell[][], ArrayList<ArrayList<Cell>>> result = iterate(grid, points_list, paths, distances_for_paths, NR_OF_ITERATIONS,
                BASE_HEIGHT_TYPE, base_scale, true, false);

        grid = result.first;
        paths = result.second;


        draw(grid, paths, NR_OF_ITERATIONS, false);

        generate_gif();

        write_output_configuration();

    }

    public static void write_output_configuration () throws IOException {
        FileWriter fileWriter = new FileWriter(currentWorkingPath.concat("\\" + storage_location_name) +"\\config.txt");
        fileWriter.write("NR_OF_ITERATIONS = " + NR_OF_ITERATIONS + "\n");
        fileWriter.write("BASE_HEIGHT_TYPE = " + BASE_HEIGHT_TYPE + "\n");
        fileWriter.write("base_scale = " + base_scale + "\n");
        fileWriter.write("height_function_scale = " + height_function_scale + "\n");
        fileWriter.write("height_function_width = " + height_function_width + "\n");
        fileWriter.close();
    }

    public static void generate_gif() throws IOException {


        File dir = new File(currentWorkingPath.concat("\\" + storage_location_name + "\\"));
        File[] files = dir.listFiles((dir1, name) -> name.endsWith(".png"));


        BufferedImage first = ImageIO.read(new File(currentWorkingPath.concat("\\" + storage_location_name + "\\image_0.png")));
        ImageOutputStream output = new FileImageOutputStream(new File(currentWorkingPath.concat("\\" + storage_location_name + "\\output_gif.gif")));

        GifSequenceWriter writer = new GifSequenceWriter(output, first.getType(), 250, true);
        writer.writeToSequence(first);


        Arrays.sort(files, new Comparator<File>(){
            @Override
            public int compare(File f1, File f2) {
                String s1 = f1.getName().substring(6,  f1.getName().indexOf("."));
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

    public static Tuple<Cell[][], ArrayList<ArrayList<Cell>>> iterate(Cell[][] grid, ArrayList points_list, ArrayList paths,
                                                                      ArrayList distances_for_paths, int NR_OF_ITERATIONS,
                                                                      String BASE_HEIGHT_TYPE, double scale, boolean save_outputs,
                                                                      boolean show_intermediate_results)
            throws IOException {

        adjust_height(grid, distances_for_paths);

        // Iterate a number of times:
        for (int i = 1; i < NR_OF_ITERATIONS; i++) {
            System.out.println("iteration : " + i);


            compute_flow(grid);

            paths = compute_paths(points_list, grid);

            if (save_outputs == true) {
                draw(grid, paths, i, show_intermediate_results);
            }

            distances_for_paths = compute_bfs(grid, paths);

            if (BASE_HEIGHT_TYPE.equals("SQRT")) {
                initialize_grid_height_sqrt(grid, scale);
            } else {
                initialize_grid_height(grid, scale);
            }
            adjust_height(grid, distances_for_paths);

        }

        compute_flow(grid);
        paths = compute_paths(points_list, grid);

        Tuple<Cell[][], ArrayList<ArrayList<Cell>>> tuple = new Tuple<Cell[][], ArrayList<ArrayList<Cell>>>(grid, paths);

        return tuple;

    }

    public static double gaussian(double x, double mu, double sigma) {

        return (double) (1.0 / (Math.sqrt(2.0 * Math.PI) * sigma) * Math.exp(-Math.pow((x - mu) / sigma, 2.0) / 2));

    }

    public static void adjust_height(Cell[][] grid, ArrayList distances_for_paths) {

        System.out.println("adjusting height");

        double[][] dist = (double[][]) distances_for_paths.get(0);

        for (int i = 0; i < nr_of_columns; i++) {

            for (int j = 0; j < nr_of_rows; j++) {

                Cell cell = grid[j][i];

                Iterator path_iterator = distances_for_paths.iterator();

                ArrayList distances_for_cell = new ArrayList();


                while (path_iterator.hasNext()) {

                    double[][] distances = (double[][]) path_iterator.next();

                    distances_for_cell.add(distances[cell.cell_y][cell.cell_x]);

                }


                height_function_width = 50;
                height_function_scale = 200000;

                double distance_1 = (double) distances_for_cell.get(0);
                double distance_2 = (double) distances_for_cell.get(1);

                double height = -height_function_scale * (gaussian(distance_1, 0, height_function_width) + gaussian(distance_2, 0, height_function_width));

//                if (distance_1 <= 2 || distance_2 <= 2) {
//                    System.out.println("before : " + grid[j][i].height);
//                }

                grid[j][i].height = grid[j][i].height + height;

//                if (distance_1 <= 2 || distance_2 <= 2) {
//
//                    System.out.println("after : " + grid[j][i].height);
//                }

            }
        }
    }

    public static void initialize_grid_height(Cell[][] grid, double scale) {

        for (int i = 0; i < nr_of_columns; i++) {
            for (int j = 0; j < nr_of_rows; j++) {

                grid[i][j].height = (0.05 * (Math.pow(grid[i][j].cell_x - source_x, 2) + Math.pow(grid[i][j].cell_y - source_y, 2)));

            }
        }
    }

    public static void initialize_grid_height_sqrt(Cell[][] grid, double scale) {

        for (int i = 0; i < nr_of_columns; i++) {
            for (int j = 0; j < nr_of_rows; j++) {

                grid[i][j].height = (0.05 * Math.sqrt(Math.pow(grid[i][j].cell_x - source_x, 2) + Math.pow(grid[i][j].cell_y - source_y, 2)));

            }
        }

    }

    public static void assign_zeros_to_grid(Cell[][] grid) {

        for (int i = 0; i < nr_of_columns; i++) {
            for (int j = 0; j < nr_of_rows; j++) {

                grid[i][j].height = (double) 0;

            }
        }
    }


    public static ArrayList compute_bfs(Cell[][] grid, ArrayList paths) {

        System.out.println("computing bfs");

        Iterator path_iterator = paths.iterator();


        // Direction vectors
        int dRow[] = {-1, 0, 1, 0};
        int dCol[] = {0, 1, 0, -1};

        ArrayList distances_for_paths = new ArrayList();

        while (path_iterator.hasNext()) {


            boolean[][] visited = new boolean[nr_of_columns][nr_of_rows];

            double[][] distances = new double[nr_of_columns][nr_of_rows];

            for (int i = 0; i < nr_of_columns; i++) {
                for (int j = 0; j < nr_of_rows; j++) {

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
        if (row < 0 || col < 0 || row >= nr_of_rows || col >= nr_of_columns)
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

        for (int i = 0; i < nr_of_columns; i++) {

            for (int j = 0; j < nr_of_rows; j++) {

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

    public static void draw_distances(Cell[][] grid, ArrayList paths, ArrayList distances_for_paths) throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(nr_of_rows, nr_of_columns);

        BufferedImage image = new BufferedImage(nr_of_rows, nr_of_columns,
                BufferedImage.TYPE_INT_ARGB);

        double[][] dist = (double[][]) distances_for_paths.get(1);

        double max_dist = dist[0][0];
        double min_dist = 0;

        for (int i = 0; i < nr_of_columns; i++) {

            for (int j = 0; j < nr_of_rows; j++) {

                if (dist[i][j] > max_dist) {
                    max_dist = dist[i][j];
                }
            }
        }

        for (int i = 0; i < nr_of_columns; i++) {
            for (int j = 0; j < nr_of_rows; j++) {

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


    public static void draw(Cell[][] grid, ArrayList paths, int image_index, boolean show_intermediate_results)
            throws IOException {

        jframe = new JFrame("panel");
        jframe.setSize(nr_of_rows, nr_of_columns);

        BufferedImage image = new BufferedImage(nr_of_rows, nr_of_columns,
                BufferedImage.TYPE_INT_ARGB);

        double max_height = grid[0][0].height;
        double min_hieght = grid[0][0].height;

        for (int i = 0; i < nr_of_columns; i++) {

            for (int j = 0; j < nr_of_rows; j++) {

                if (grid[i][j].height > max_height) {
                    max_height = grid[i][j].height;
                }
                if (grid[i][j].height < min_hieght) {
                    min_hieght = grid[i][j].height;
                }
            }
        }

        System.out.println("min : " + min_hieght + " max: " + max_height);

        for (int i = 0; i < nr_of_columns; i++) {
            for (int j = 0; j < nr_of_rows; j++) {

                int c = (int) (grid[i][j].height * 255.0 / max_height);


                if (c < 0) {
                    c = (int) (grid[i][j].height * 255.0 / min_hieght);
                    image.setRGB(i, j, new Color(0, 0, c).getRGB());
                } else {
                    image.setRGB(i, j, new Color(c, 0, 0).getRGB());
                }


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

        File file = new File(currentWorkingPath.concat("/" + storage_location_name + "/image_" + image_index + ".png"));
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

            if (point.name.equals(target_name)) {
                continue;
            }

            Cell grid_cell = grid[point.grid_x][point.grid_y];

            Cell current_cell = grid_cell;

            ArrayList<Cell> path = new ArrayList();

            path.add(current_cell);

            while (!(current_cell.title.equals(target_name))) {

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
        for (int i = 0; i < nr_of_columns; i++) {

            for (int j = 0; j < nr_of_rows; j++) {

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

                if (y + 1 < nr_of_rows) {
                    bottom = grid[x][y + 1];
                    neighbors.add(bottom);
                }

                if (x + 1 < nr_of_columns) {
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

                if (y - 1 >= 0 && x + 1 < nr_of_columns) {
                    top_right = grid[x + 1][y - 1];
                    neighbors.add(top_right);
                }

                if (x - 1 >= 0 && y + 1 < nr_of_rows) {
                    bottom_left = grid[x - 1][y + 1];
                    neighbors.add(bottom_left);
                }

                if (x + 1 < nr_of_columns && y + 1 < nr_of_rows) {
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

            if (point.name.equals(target_name)) {
                grid[point.grid_x][point.grid_y].height = 0;

                source_x = point.grid_x;
                source_y = point.grid_y;

            }
        }
    }


    public static void compute_cell_for_point(Bounds bounds, ArrayList points) {

        Iterator it = points.iterator();

        double extent_x = bounds.max_x - bounds.min_x;
        double extent_y = bounds.max_y - bounds.min_y;

        double step_x = extent_x / nr_of_columns;
        double step_y = extent_y / nr_of_rows;

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

        Scanner sc = new Scanner(new File("./input/1_s_2_t.csv"));
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
