from copy import deepcopy
import sympy as sp
from sympy import Symbol


# returns the m (slope) between 2 points
def find_slope(p1: tuple, p2: tuple):
    y = p2[1] - p1[1]
    x = p2[0] - p1[0]
    return y / x


# returns an identity matrix based on n size
def identity_matrix(n):
    matrix = []
    for i in range(n):
        matrix.append([0 for _ in range(n)])
    for j in range(n):
        matrix[j][j] = 1
    return matrix


# returns the inverted A matrix
def invert_matrix(A):
    A_COPY = deepcopy(A)
    I = identity_matrix(len(A))
    indices = list(range(len(A)))

    for fd in range(len(A)):
        pivot = 1 / A_COPY[fd][fd]
        for j in range(len(A)):
            A_COPY[fd][j] *= pivot
            I[fd][j] *= pivot

        for i in indices[0:fd] + indices[fd + 1:]:
            crScaler = A_COPY[i][fd]
            for j in range(len(A)):
                A_COPY[i][j] = A_COPY[i][j] - crScaler * A_COPY[fd][j]
                I[i][j] = I[i][j] - crScaler * I[fd][j]
    return I


# returns the result of polynom for n,m based on the formula
def get_polynom_result(points: list, m: int, n: int):
    if m + 1 == n:
        return points[n][1]
    if m == n - 1:
        return points[m][1]


# returns the L(x) for position i
def find_L(points: list, i: int, x: Symbol):
    f = 1
    for j in range(len(points)):
        # only if i != j we add the next equation
        if j != i:
            f *= (x - points[j][0]) / (points[i][0] - points[j][0])
    return f


# returns the solution that derives from multiplying matrix * vector
# we know that (1, x, x^2) vector has the length of 3,
# so we multiply each item in the row with the matching item in the vector
def solve_matrix(matrix: list, vector: list):
    solutions = []
    s = 0
    for row in matrix:
        for i in range(3):
            s += (vector[i] * row[i])
        solutions.append(s)
        s = 0
    return solutions


# returns a lambdified polynom from the list of values
# iterates over each value and multiply it by the current degree
# for example: value #0 => value * x ** 0
def create_polynom(values: list, x: Symbol):
    f = 0
    degree = 0
    for value in values:
        f += (value * x ** degree)
        degree += 1

    return sp.lambdify(x, f)


# returns the result using the linear method
# takes 2 points which our value is in between, and creates a linear equation
# from that equation we return the result by assigning the value
def linear_method(points: list, value):
    result = None
    for i in range(len(points)):
        if points[i][0] > value:
            if i - 1 < 0:
                print('\033[91m' + "Cannot find a linear equation! "
                                   "value is not between 2 valid points.\nExiting..")
                exit(1)

            point1 = points[i - 1]
            point2 = points[i]
            m = find_slope(point1, point2)
            x = sp.symbols('x')
            f = m * x + (point1[1] - m * point1[0])
            result = sp.lambdify(x, f)(value)

    return result


# returns the result using the polynomial method
# creates a polynom using the points list
# from the polynom created, we assign the value and return the result
def polynomial_method(points: list, value):
    solutions = [point[1] for point in points]  # vector of Y's
    matrix = [[1, point[0], point[0] ** 2] for point in points]

    if len(matrix) != len(matrix[0]):
        print('\033[91m' +
              "Cannot create an inverted matrix to non square matrix!\nExiting..")
        exit(1)

    values = solve_matrix(invert_matrix(matrix), solutions)

    # after we found the values, now we create the polynom
    x = sp.symbols('x')
    polynom = create_polynom(values, x)
    return polynom(value)


# returns the result using the lagrange method
# we generate a polynom created from the sum of Li(x) using the formula
# after the polynom is created, we assign the value and return the result
def lagrange_method(points: list, value):
    f = 0
    x = sp.symbols('x')
    for i in range(len(points)):
        l = find_L(points, i, x)
        f += l * points[i][1]  # Li(x) * Yi

    return sp.lambdify(x, f)(value)


# returns the result using the lagrange method
# we find the result of each polynom created from the list 1..n
# we use our previous solutions to find the largest polynom p1->n
def neville_method(points: list, value):
    solutions = {}
    x = sp.symbols('x')
    size = len(points)

    # first, we add the consecutive pairs such as 01, 12, 23 ..
    # after we add the solutions to the dictionary, we can use
    # the results we found to proceed to 02, 13 etc..
    for i in range(size - 1):
        m = i
        n = i + 1
        x1 = points[m][0]
        x2 = points[n][0]
        f = ((x - x1) * points[n][1] - (x - x2) * points[m][1]) / (x2 - x1)
        solutions[f"{m}{n}"] = sp.lambdify(x, f)(value)

    # once we have to consecutive pairs, we can continue based on the current interval
    # we move forward until the inteval equals to the size of the list
    # e.g: interval = 2 => 02 13 24 ...  interval = 3 => 03 14 25 ...
    interval = 2
    for i in range(size):
        if interval == size:
            break

        for j in range(size - interval):
            m = j
            n = j + interval
            x1 = points[m][0]
            x2 = points[n][0]
            f = ((x - x1) * solutions[f"{m + 1}{n}"] - (x - x2)
                 * solutions[f"{m}{n - 1}"]) / (x2 - x1)
            solutions[f"{m}{n}"] = sp.lambdify(x, f)(value)
        interval += 1

    # return the result of the polynom from the beginning index to the last index
    return solutions[f"0{size - 1}"]


# the main method
def main():
    selected_method = None
    points = [(1, 0.7651), (1.3, 0.62), (1.6, 0.4554)]
    x = 1.5  # the x we will recieve the f(x) to

    print("Choose method: \n1) Linear\n2)Polynomial\n3)Lagrange\n4)Neville")
    option = input("Choose option: ")

    # input check for invalid option
    while not option.isdigit() or (int(option) < 1 or int(option) > 4):
        option = input("Incorrect option! Try again!\nChoose option: ")
    option = int(option)

    if option == 1:
        selected_method = linear_method
    if option == 2:
        selected_method = polynomial_method
    if option == 3:
        selected_method = lagrange_method
    if option == 4:
        selected_method = neville_method

    print(f"Selected method: {selected_method.__name__.replace('_', ' ')}")
    result = selected_method(points, x)
    print(f"the result is: {result}")


if __name__ == '__main__':
    main()
