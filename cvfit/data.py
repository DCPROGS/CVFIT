

class DataPoint(object):
    """
    """

    def __init__(self, x, y, se=0):
        self.x = x
        self.y = y
        self.se = se
        if se > 0:
            self.w = 1.0 / (se * se)
        else:
            self.w = 1.0

class DataSet(object):
    """
    """
    def __init__(self):
        self.points = []

    def from_columns(self, x, y):
        for i in range(len(x)):
            self.points.append(DataPoint(x[i], y[i]))

    def add_point(self, x, y):
        self.points.append(DataPoint(x, y))

    def size(self):
        return len(self.points)

    def __str__(self):
        str = "\n"
        for point in self.points:
            str += "{0:.6g}\t{1:.6g}\n".format(point.x, point.y)
        return str