import tkinter as tk
from tkinter import messagebox
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import schemdraw
import schemdraw.elements as elm
from matplotlib.animation import FuncAnimation

class RLC_Circuit:
    def __init__(self, R, L, C, V_func):
        self.R = R
        self.L = L
        self.C = C
        self.V_func = V_func

    def equations(self, t, y):
        q, i = y
        dqdt = i
        didt = (self.V_func(t) - self.R * i) / self.L
        return [dqdt, didt]

    def simulate(self, t_span, y0, t_eval=None):
        if t_eval is None:
            t_eval = np.linspace(*t_span, 1000)
        solution = solve_ivp(
            self.equations, t_span, y0, t_eval=t_eval, method="RK45"
        )
        return {
            "time": solution.t,
            "charge": solution.y[0],
            "current": solution.y[1],
            "voltage_capacitor": solution.y[0] / self.C,  # Correct voltage across the capacitor
        }

    def draw_circuit(self):
        """
        Draw a schematic of the RLC circuit with the capacitor connected in parallel to the voltage source,
        and the resistor and inductor connected in series.
        """
        with schemdraw.Drawing() as d:
            # Drawing components

            V1 = elm.SourceV().up().label(f'V = {self.V_func(0)}', loc='top')
            elm.Switch().right()
            R = elm.Resistor().label(f'R = {self.R} Ω', loc='top')
            Q = elm.Diode().down()
            I = elm.Inductor().left().label(f'L = {self.L} H', loc='bottom')
            elm.Ground()

            elm.Line().tox(V1.start)



            elm.Capacitor().endpoints(I.end, R.start).label(f'C = {self.C} F', loc='bottom')

            d.draw()  # Draw the schematic

# GUI for user inputs
def run_gui():
    def simulate_circuit():
        try:
            print("Starting simulation...")  # Debugging line
            R = float(resistance_entry.get())
            L = float(inductance_entry.get())
            C = float(capacitance_entry.get())
            V_value = float(voltage_entry.get())

            # Define the voltage function
            V_func = lambda t: V_value

            # Create an RLC circuit object
            rlc_circuit = RLC_Circuit(R, L, C, V_func)

            # Now include the transient current animation as part of the simulation
            animate_transient_response(R, L, C, V_value)

        except ValueError as e:
            print(f"Error: {e}")  # Print error for debugging
            messagebox.showerror("Invalid Input", "Please enter valid numerical values.")

    def generate_schematic():
        try:
            R = float(resistance_entry.get())
            L = float(inductance_entry.get())
            C = float(capacitance_entry.get())
            V_value = float(voltage_entry.get())

            # Define the voltage function
            V_func = lambda t: V_value

            # Create an RLC circuit object
            rlc_circuit = RLC_Circuit(R, L, C, V_func)

            # Draw the circuit schematic
            rlc_circuit.draw_circuit()

        except ValueError as e:
            print(f"Error: {e}")  # Print error for debugging
            messagebox.showerror("Invalid Input", "Please enter valid numerical values.")

    def animate_transient_response(R, L, C, V):
        # Define the transient time range
        t_values = np.linspace(0, 2, 500)  # Time range limited to show transient response

        if R** - (4 * L) / C <= 0:
            print("Not Overdamped")

        else:
            a = L * C
            b = R * C
            c = 1

            m1 = (-b + (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
            m2 = (-b - (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)

            c1 = m1
            c2 = m2

            l1 = c1-c2

            k1 = V/l1
            k2 = -k1


            # Define the components of the response
            i1_values = k1 * np.exp(m1 * t_values)  # First exponential component
            i2_values = k2 * np.exp(m2 * t_values)  # Second exponential component
            i_total = i1_values + i2_values  # Total response

            # Set up the figure and axis for animation
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_xlim(0, 2)
            ax.set_ylim(min(i_total) - 0.25, max(i_total) + 0.25)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("i(t) (amperes)")
            ax.set_title("Transient Current of an Overdamped RLC Circuit")
            ax.grid(True, linestyle='--')

            # Plot the static components
            ax.plot(t_values, i1_values, label=r"$+1.25 e^{-t}$", linestyle='--', color='gray')
            ax.plot(t_values, i2_values, label=r"$-1.25 e^{-9t}$", linestyle='--', color='black')

            # Initialize the total response line
            line, = ax.plot([], [], label="Total Response", color='teal', linewidth=2)
            ax.legend()

        # Initialize the animation
        def init():
            line.set_data([], [])
            return line,

        # Update function for the animation
        def update(frame):
            line.set_data(t_values[:frame], i_total[:frame])
            return line,

        # Create the animation
        ani = FuncAnimation(fig, update, frames=len(t_values), init_func=init, blit=True, interval=30)

        # Show the plot with animation
        plt.show()

    # GUI Window
    window = tk.Tk()
    window.title("RLC Circuit Simulator")

    # Input fields
    tk.Label(window, text="Resistance (Ω):").grid(row=0, column=0)
    resistance_entry = tk.Entry(window)
    resistance_entry.grid(row=0, column=1)

    tk.Label(window, text="Inductance (H):").grid(row=1, column=0)
    inductance_entry = tk.Entry(window)
    inductance_entry.grid(row=1, column=1)

    tk.Label(window, text="Capacitance (F):").grid(row=2, column=0)
    capacitance_entry = tk.Entry(window)
    capacitance_entry.grid(row=2, column=1)

    tk.Label(window, text="Voltage Source (V):").grid(row=3, column=0)
    voltage_entry = tk.Entry(window)
    voltage_entry.grid(row=3, column=1)

    # Simulate button
    simulate_button = tk.Button(window, text="Simulate", command=simulate_circuit)
    simulate_button.grid(row=4, columnspan=2)

    # Generate schematic button
    schematic_button = tk.Button(window, text="Generate Schematic", command=generate_schematic)
    schematic_button.grid(row=5, columnspan=2)

    window.mainloop()

# Entry point of the script
if __name__ == "__main__":
    run_gui()
