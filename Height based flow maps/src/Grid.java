import java.util.ArrayList;
import java.util.Iterator;

public class Grid {

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
    public static Bounds obtain_bounds(ArrayList<Point> input_points) {

        Iterator it = input_points.iterator();

        double min_x = input_points.get(0).col;
        double min_y = input_points.get(0).row;
        double max_x = input_points.get(0).col;
        double max_y = input_points.get(0).row;

        while (it.hasNext()) {

            Point next_point = (Point) it.next();

            if (next_point.col < min_x) {
                min_x = next_point.col;
            }
            if (next_point.col > max_x) {
                max_x = next_point.col;
            }
            if (next_point.row < min_y) {
                min_y = next_point.row;
            }
            if (next_point.row > max_y) {
                max_y = next_point.row;
            }

        }

        Bounds result = new Bounds();
        result.max_x = max_x;
        result.min_x = min_x;
        result.min_y = min_y;
        result.max_y = max_y;

        return result;
    }

    public static void initialize_points_in_grid(Cell[][] grid, ArrayList<Point> points_list, String TARGET_NAME, int source_x, int source_y, Cell source_cell) {

        Iterator it = points_list.iterator();

        // setting the names of source/targets and setting height of source
        while (it.hasNext()) {

            Point point = (Point) it.next();

            grid[point.grid_col][point.grid_row].title = point.name;

            if (point.name.equals(TARGET_NAME)) {
                grid[point.grid_col][point.grid_row].height = -10000;

                source_x = point.grid_col;
                source_y = point.grid_row;

            }
        }

        source_cell = grid[source_x][source_y];
    }

}
