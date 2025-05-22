# Python Basics for OptVL Users

## Introduction

Welcome to OptVL! If you're here, you're likely interested in leveraging the power of Athena Vortex Lattice (AVL) within the flexible Python environment. You might be wondering if you need to be a Python expert to use OptVL. The short answer is: **no, you don't!**

This guide is designed to give you the bare minimum Python knowledge needed to understand and run the OptVL examples. OptVL itself is a Python library, which means all the examples and scripts you'll encounter are written in Python. Our goal here is to make you comfortable with the fundamental concepts you'll see in these examples.

## 1. What is a Python Script?

A Python script is simply a plain text file that contains a sequence of instructions written in the Python language. These files typically have a `.py` extension. For example, you might see scripts named `run_my_analysis.py` or `plot_results.py`.

When you "run" a Python script, you're telling the Python interpreter (the program that understands Python) to read the instructions in your file and execute them one by one.

## 2. How to Run OptVL Examples/Python Scripts

Before you can run any Python scripts, including the OptVL examples, you need to have Python and OptVL installed. If you haven't done this yet, please refer to the OptVL installation guide.

Once everything is set up, the most common way to run a Python script is from your command line or terminal. You typically navigate to the directory where the script is located and then type:

```bash
python your_script_name.py
```

For example, to run an OptVL example named `run_basic.py`, you would type:

```bash
python run_basic.py
```

**A Note on IDEs:**
Some users prefer to work with Integrated Development Environments (IDEs) like Spyder, VS Code, or PyCharm. These tools offer more features, such as code editors, debuggers, and often a "run" button that executes the current script. While IDEs can be very helpful, this guide will focus on the command-line approach for its simplicity and universality. The core Python concepts remain the same regardless of how you run the script.

## 3. Basic Python Syntax and Concepts in OptVL Examples

Let's dive into the essential Python elements you'll encounter when working with OptVL examples.

### Comments

In Python, the `#` symbol is used to denote a comment. Anything on a line after a `#` is ignored by the Python interpreter. Comments are incredibly useful for explaining what the code is doing, both for yourself and for others who might read your code.

```python
# This is a comment. It explains the next line of code.
alpha = 5.0 # This sets the angle of attack to 5 degrees.
```

### Variables

Variables are names that you give to store data in your program. Think of them as labels for values. You assign a value to a variable using the equals sign (`=`).

```python
# Examples of variables
angle_of_attack = 5.0
aircraft_filename = "my_plane.avl"
is_analysis_complete = True # A boolean variable (True/False)
```

In OptVL examples, you'll commonly see variables holding different types of data:

*   **Strings:** These are sequences of text characters, enclosed in either single quotes (`'`) or double quotes (`"`).
    ```python
    wing_name = "main_wing"
    input_file = 'aircraft_geometry.avl'
    ```

*   **Numbers:** These can be integers (whole numbers) or floating-point numbers (numbers with a decimal point).
    ```python
    number_of_cases = 10      # An integer
    mach_number = 0.75        # A floating-point number
    reference_area = 30.5     # Another float
    ```

*   **Lists:** These are ordered collections of items, enclosed in square brackets `[]`, with items separated by commas. Lists are useful for storing multiple related values.
    ```python
    angles = [0.0, 2.0, 4.0, 6.0]
    control_surfaces = ["elevator", "aileron", "rudder"]
    ```

*   **Dictionaries:** These are unordered collections of key-value pairs, enclosed in curly braces `{}`. Each item consists of a key (often a string) and its associated value, separated by a colon `:`. Dictionaries are very common for accessing results from OptVL.
    ```python
    # Example of a dictionary
    parameters = {'mach': 0.8, 'altitude': 10000}
    results = {'CL': 0.525, 'CD': 0.051, 'Cm': -0.02}

    # You access values using their keys:
    lift_coefficient = results['CL'] # This would be 0.525
    drag_coefficient = results['CD']   # This would be 0.051
    print(lift_coefficient)
    ```

### Importing Libraries

Python has a vast collection of libraries that provide additional functionality. OptVL is one such library. To use a library in your script, you first need to `import` it.

You'll typically see this at the beginning of OptVL example scripts:

```python
import optvl
```

Sometimes, you might see a variation like:

```python
from optvl import OptVL # Imports only the OptVL class
```

Either way, the `import` statement makes the functions, classes, and other tools provided by OptVL available for you to use in your script.

### Creating an OptVL Object

After importing OptVL, you'll usually create an "OptVL object." This object will represent your AVL session and is what you'll use to interact with AVL's functionalities.

A common way to do this is:

```python
import optvl

# Create an OptVL object
avl_session = optvl.OptVL()
```
Or, if you used `from optvl import OptVL`:
```python
from optvl import OptVL

# Create an OptVL object
avl_session = OptVL()
```

Now, the variable `avl_session` (you can name it anything, but `avl` or `avl_session` is common) holds your OptVL instance.

### Function/Method Calls

Much of your interaction with OptVL will be through calling "functions" or "methods." Methods are functions that belong to an object. The syntax is typically `object_name.method_name(argument1, argument2, ...)`.

Arguments are the inputs you provide to the method to tell it what to do.

```python
# Example: Load an AVL geometry file
avl_session.load_geometry("aircraft.avl")

# Example: Set an aerodynamic parameter (angle of attack to 5 degrees)
avl_session.set_parameter("alpha", 5.0)

# Example: Execute the analysis run
results = avl_session.execute_run()
```
In these examples:
*   `load_geometry`, `set_parameter`, and `execute_run` are methods of the `avl_session` object.
*   `"aircraft.avl"`, `"alpha"`, and `5.0` are arguments passed to these methods.
*   `execute_run()` returns some data (the analysis results), which we store in a variable named `results`.

### Printing Output

The `print()` function is a built-in Python function used to display information to your command line or terminal. This is very useful for seeing results, status messages, or debugging.

```python
print("Starting the aerodynamic analysis...")
CL = 0.65
print("Calculated Lift Coefficient (CL):", CL) # You can print variables
print(results) # You can print complex objects like dictionaries
```

### Simple Loops

You might encounter `for` loops in some OptVL examples, especially when performing tasks for multiple values (e.g., running an analysis for a list of angles of attack). A `for` loop allows you to repeat a block of code for each item in a sequence (like a list).

Here's a basic structure:

```python
# A list of angles of attack
alpha_values = [0.0, 2.5, 5.0, 7.5]

# Loop through each angle in the list
for current_alpha in alpha_values:
    # The code inside this loop (indented) will run for each angle
    avl_session.set_parameter("alpha", current_alpha)
    current_results = avl_session.execute_run()
    print("Results for alpha =", current_alpha, ":")
    print("CL =", current_results['CL'], "CD =", current_results['CD'])

print("Loop finished.")
```
In this example:
*   The indented block of code (setting alpha, executing run, printing results) is executed once for each value in `alpha_values`.
*   In the first iteration, `current_alpha` will be `0.0`. In the second, it will be `2.5`, and so on.

Indentation (the spaces at the beginning of lines) is very important in Python. It's how Python knows which lines of code are inside the loop and which are outside.

## 4. Understanding OptVL Examples

The best way to get comfortable is to open up the OptVL example scripts (like `examples/run_basic.py` or `examples/run_aero_sweeps.py`). Try to identify the elements we've discussed:
*   Comments explaining the code.
*   Variables storing filenames, parameters, or results.
*   `import optvl` at the beginning.
*   The creation of an OptVL object (e.g., `avl = optvl.OptVL()`).
*   Method calls on the OptVL object (e.g., `avl.load_geometry(...)`, `avl.set_parameter(...)`, `avl.execute_run()`).
*   `print()` statements showing output.
*   Any `for` loops used to iterate through cases.

**Try This:**
Once you've identified these parts, try making small modifications to an example script. For instance:
*   Change an angle of attack value.
*   Change the name of the AVL input file (if you have another one).
*   Add more `print()` statements to see the values of different variables.
Then, run the script to see how your changes affect the outcome. This hands-on approach is a great way to learn!

## 5. Where to Learn More Python

As we said at the start, you don't need to become a Python guru to use OptVL effectively. The concepts covered here should be sufficient for understanding and running the provided examples and making basic modifications.

However, if you find yourself wanting to write more complex scripts, customize OptVL's behavior significantly, or simply become more proficient in Python, here are a couple of excellent (and free) resources for beginners:

*   **The Official Python Tutorial:** [https://docs.python.org/3/tutorial/](https://docs.python.org/3/tutorial/)
    *   Comprehensive and authoritative, directly from the creators of Python.
*   **W3Schools Python Tutorial:** [https://www.w3schools.com/python/](https://www.w3schools.com/python/)
    *   Known for its clear explanations and interactive "try it yourself" examples.

Many other fantastic tutorials and courses are available online if these don't suit your learning style.

## Conclusion

With these Python basics under your belt, you should be well-equipped to dive into the OptVL examples and start exploring its capabilities. Remember that the examples are there to be experimented with – don't be afraid to change things and see what happens!

Now, you should be ready to proceed to tutorials that focus on specific OptVL functionalities, knowing the basic Python syntax that underpins them. Happy analyzing!
