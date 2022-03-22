import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class Main {


    public static int nr_of_rows = 100;
    public static int nr_of_columns = 100;

    public static int source_x = 0;
    public static int source_y = 0;

    public static String target_name = "A";

    public static void main(String[] args) throws FileNotFoundException {

        System.out.println("test");

        //Cell cell = new Cell();

        //int nr_of_rows = 100;
        //int nr_of_columns = 100;

        //String target_name = "A";

        Cell[][] grid = initialize_grid(nr_of_rows, nr_of_columns);


        ArrayList items = read_input();

        ArrayList points_list = process_input(items);

        Bounds bounds = obtain_bounds(points_list);

        compute_cell_for_point(bounds, points_list);


        assign_height_to_grid(grid);


        System.out.println();


    }

    public static void initialize_points_in_grid(Cell[][] grid, ArrayList<Point> points_list) {

        Iterator it = points_list.iterator();

        // setting the names of source/targets and setting height of source
        while (it.hasNext()) {

            Point point = (Point) it.next();

            grid[point.grid_x][point.grid_y].title = point.name;

            if (point.name.equals(target_name)) {
                grid[point.grid_x][point.grid_y].height = -20;

                source_x = point.grid_x;
                source_y = point.grid_y;

            }
        }
    }


    public static void assign_height_to_grid(Cell[][] grid) {

        for (int i = 0; i < nr_of_columns; i++) {
            for (int j = 0; j < nr_of_rows; j++) {

                grid[i][j].height = (float) (grid[i][j].height + 0.05 * (Math.pow(grid[i][j].cell_x - source_x, 2) + Math.pow(grid[i][j].cell_y - source_y, 2)));

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

        float min_x = input_points.get(0).x;
        float min_y = input_points.get(0).y;
        float max_x = input_points.get(0).x;
        float max_y = input_points.get(0).y;

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

        //parsing a CSV file into Scanner class constructor
        //System.out.println("Working Directory = " + System.getProperty("user.dir"));

        Scanner sc = new Scanner(new File("./input/1_s_2_t.csv"));
        sc.useDelimiter("\n");   //sets the delimiter pattern

        ArrayList items = new ArrayList();

        while (sc.hasNext()) {
            //System.out.print(sc.next());

            items.add(sc.next());

        }
        sc.close();  //closes the scanner

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

                cell.title = Integer.toString(i) + " " + Integer.toString(j);

                grid[i][j] = cell;

                //System.out.println(grid[i][j].temp);
                //grid[i][j]
            }
        }

        return grid;
    }

}
