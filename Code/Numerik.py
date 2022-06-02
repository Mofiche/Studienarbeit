from sklearn.metrics import r2_score
import numpy as np


class RungeKutta4(object):

    GENAUIGKEIT = 1000

    def DGL_Test(self, x, y):
        dy = y
        return dy

    def solve_DGL_RK4(self, DGL, xMIN, xMAX, y0):
        xi = xMIN
        yi = y0
        x, y = [], []
        diff = []
        h = (xMAX - xMIN) / self.GENAUIGKEIT
        print("xMIN", xMIN, "yMin", y0)

        for _ in range(self.GENAUIGKEIT):
            x.append(xi)
            y.append(yi)
            k1 = DGL(xi, yi)
            k2 = DGL(xi + (h / 2), yi + (h / 2) * k1)
            k3 = DGL(xi + (h / 2), yi + (h / 2) * k2)
            k4 = DGL(xi + h, yi + k3)
            dy = (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            diff.append(dy)
            yi += dy
            xi += h

        print("xMAX", xi, "yMax", yi)
        """
        coeffs_calc = np.polyfit(x, y, 50)
        mymodel_calc = np.poly1d(coeffs_calc)

        coeffs = np.polyfit(x, y, 5)
        mymodel = np.poly1d(coeffs)

        R = "R = {:.2f} %".format(100 * r2_score(y, mymodel(x)))
        R_calc = "R_calc = {:.2f} %".format(100 * r2_score(y, mymodel_calc(x)))
        """
        return x, y, diff

    def solve_DGL_RK4_InputArray(self, DGL, x_input, y0):
        xi = min(x_input)
        xmin = xi
        xmax = max(x_input)
        yi = y0
        x, y = [], []
        h = (xmax - xmin) /len(x_input)
        print("xMIN", xmin, "yMin", y0)

        for _ in range(len(x_input)):
            x.append(xi)
            y.append(yi)
            k1 = DGL(xi, yi)
            k2 = DGL(xi + (h / 2), yi + (h / 2) * k1)
            k3 = DGL(xi + (h / 2), yi + (h / 2) * k2)
            k4 = DGL(xi + h, yi + k3)

            yi += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            xi += h

        print("xMAX", xi, "yMax", yi)

        return y