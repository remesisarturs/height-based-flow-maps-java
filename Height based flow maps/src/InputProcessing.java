import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class InputProcessing {

    public static ArrayList read_input(String INPUT_FILE_NAME) throws FileNotFoundException {

        Scanner sc = new Scanner(new File(INPUT_FILE_NAME));
        sc.useDelimiter("\n");

        ArrayList items = new ArrayList();

        while (sc.hasNext()) {

            items.add(sc.next());

        }
        sc.close();

        return items;
    }

    public static ArrayList process_input(ArrayList items) {

        Iterator it = items.iterator();

        ArrayList point_list = new ArrayList();

        while (it.hasNext()) {

            String input_line = it.next().toString();
            List<String> components = Arrays.asList(input_line.split(";"));

            Point point = new Point();
            point.name = components.get(0);
            point.col = Float.parseFloat(components.get(1));
            point.row = Float.parseFloat(components.get(2));

            point_list.add(point);

        }

        return point_list;
    }

    public static void compute_cell_for_point(Bounds bounds, ArrayList points, int NR_OF_COLUMNS, int NR_OF_ROWS) {

        Iterator it = points.iterator();

        double extent_x = bounds.max_x - bounds.min_x;
        double extent_y = bounds.max_y - bounds.min_y;

        double step_x = extent_x / NR_OF_COLUMNS;
        double step_y = extent_y / NR_OF_ROWS;

        while (it.hasNext()) {

            Point point = (Point) it.next();

            int i = 1;
            int j = 1;

            while (bounds.min_x + i * step_x < point.col) {
                i = i + 1;
            }
            int column_index = i - 1;

            while (bounds.min_y + j * step_y < point.row) {
                j = j + 1;
            }
            int row_index = j - 1;

            point.grid_col = column_index;
            point.grid_row = row_index;

        }
    }

}
