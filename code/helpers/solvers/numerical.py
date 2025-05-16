from sympy import Expr, Symbol, diff, lambdify
from typing import Union

def newtons_method(function:Expr, variable:Symbol, guess: Union[float, int], tolerance=1e-7, max_iterations=1000) -> float:
    """
    Newton's Method to numerically solve for the root of an algebraic function.

    Parameters:
    function (symbolic): The function whose root is to be found.
    variable (symbol): The independent variable of the function.
    guess (float/int): Initial guess for the root.
    tolerance (float): The accuracy of the solution.
    max_iter (int): Maximum number of iterations.
    """
    x = variable
    f = function
    df = diff(f, x)
    f_, df_ = lambdify(x, f), lambdify(x, df)
    x_now = guess
    for i in range(max_iterations):
        fx = f_(x_now)
        dfx = df_(x_now)

        if abs(dfx) < 1e-10: # Avoid division by zero
            raise ValueError("Derivative near zero; Newton's method fails.")
        
        x_new = x_now - fx / dfx

        print(f"Iteration {i}: x = {x_now}, f(x) = {fx}, f'(x) = {dfx}")

        if abs(x_new - x_now) < tolerance:
            return x_new
        x_now = x_new

    raise ValueError("Newton's method did not converge within the maximum number of iterations")
