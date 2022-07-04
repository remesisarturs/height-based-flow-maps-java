import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class InputProcessing {

    public static ArrayList readInput(String INPUT_FILE_NAME) throws FileNotFoundException {

        Scanner sc = new Scanner(new File(INPUT_FILE_NAME));
        sc.useDelimiter("\n");

        ArrayList items = new ArrayList();

        while (sc.hasNext()) {

            items.add(sc.next());

        }
        sc.close();

        return items;
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
            point.row = Float.parseFloat(components.get(2));

            pointList.add(point);

        }

        return pointList;
    }

    public static void computeCellForPoint(Bounds bounds, ArrayList points, int NR_OF_COLUMNS, int NR_OF_ROWS) {

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

}
