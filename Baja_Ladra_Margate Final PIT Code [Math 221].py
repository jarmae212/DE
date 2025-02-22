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
            "voltage_capacitor": solution.y[0] / self.C,  
        }

    def draw_circuit(self):
        """
        Draw a schematic of the RLC circuit with the capacitor connected in parallel to the voltage source,
        and the resistor and inductor connected in series.
        """
        with schemdraw.Drawing() as d:
            
            V1 = elm.SourceV().up().label(f'V = {self.V_func(0)}', loc='top')
            elm.Switch().right()
            R = elm.Resistor().label(f'R = {self.R} Ω', loc='top')
            Q = elm.Diode().down()
            I = elm.Inductor().left().label(f'L = {self.L} H', loc='bottom')
            elm.Ground()

            elm.Line().tox(V1.start)
            elm.Capacitor().endpoints(I.end, R.start).label(f'C = {self.C} F', loc='bottom')

            d.draw()  

def run_gui():
    def simulate_circuit():
        try:
            print("Starting simulation...")  
            R = float(resistance_entry.get())
            L = float(inductance_entry.get())
            C = float(capacitance_entry.get())
            V_value = float(voltage_entry.get())

           
            V_func = lambda t: V_value
       
            rlc_circuit = RLC_Circuit(R, L, C, V_func)

            animate_transient_response(R, L, C, V_value)

        except ValueError as e:
            print(f"Error: {e}") 
            messagebox.showerror("Invalid Input", "Please enter valid numerical values.")

    def generate_schematic():
        try:
            R = float(resistance_entry.get())
            L = float(inductance_entry.get())
            C = float(capacitance_entry.get())
            V_value = float(voltage_entry.get())

            
            V_func = lambda t: V_value

            
            rlc_circuit = RLC_Circuit(R, L, C, V_func)

            
            rlc_circuit.draw_circuit()

        except ValueError as e:
            print(f"Error: {e}")   
            messagebox.showerror("Invalid Input", "Please enter valid numerical values.")

    def animate_transient_response(R, L, C, V):
        if R**2 == (4 * L )/ C:
            print("Critically Damped Circuit")

            alpha = R / (2 * L)  
            A = V / L  

            t_values = np.linspace(0, 2, 500)

            i_values = A * t_values * np.exp(-alpha * t_values)
            
            i_upper_values = A * t_values * np.exp(-alpha * t_values)
            i_lower_values = -A * t_values * np.exp(-alpha * t_values)

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_xlim(0, 2)
            ax.set_ylim(min(i_lower_values) - 0.25, max(i_upper_values) + 0.25)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("i(t) (amperes)")
            ax.set_title("Critically Damped Current of an RLC Circuit")
            ax.grid(True, linestyle='--')

            line, = ax.plot([], [], label="Total Response", color='blue', linewidth=3)
            
            upper_envelope, = ax.plot([], [], label="Upper Envelope", color='orange', linestyle='--', linewidth=2)
            lower_envelope, = ax.plot([], [], label="Lower Envelope", color='red', linestyle='--', linewidth=2)
            
            ax.legend()

            def init():
                line.set_data([], [])
                upper_envelope.set_data([], [])
                lower_envelope.set_data([], [])
                return line, upper_envelope, lower_envelope

            def update(frame):
                line.set_data(t_values[:frame], i_values[:frame])
                upper_envelope.set_data(t_values[:frame], i_upper_values[:frame])
                lower_envelope.set_data(t_values[:frame], i_lower_values[:frame])
                return line, upper_envelope, lower_envelope

            ani = FuncAnimation(fig, update, frames=len(t_values), init_func=init, blit=True, interval=30)

            plt.show()

        else:
            messagebox.showerror("Not Critically Damped", "This is not a critically damped circuit. Check R, L, and C values.")

    window = tk.Tk()
    window.title("RLC Circuit Simulator")

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

    simulate_button = tk.Button(window, text="Simulate", command=simulate_circuit)
    simulate_button.grid(row=4, columnspan=2)

    schematic_button = tk.Button(window, text="Generate Schematic", command=generate_schematic)
    schematic_button.grid(row=5, columnspan=2)

    window.mainloop()

if __name__ == "__main__":
    run_gui()
    inductance_entry.grid(row=1, column=1)

    tk.Label(window, text="Capacitance (F):").grid(row=2, column=0)
    capacitance_entry = tk.Entry(window)
    capacitance_entry.grid(row=2, column=1)

    tk.Label(window, text="Voltage Source (V):").grid(row=3, column=0)
    voltage_entry = tk.Entry(window)
    voltage_entry.grid(row=3, column=1)

    simulate_button = tk.Button(window, text="Simulate", command=simulate_circuit)
    simulate_button.grid(row=4, columnspan=2)

    schematic_button = tk.Button(window, text="Generate Schematic", command=generate_schematic)
    schematic_button.grid(row=5, columnspan=2)

    window.mainloop()

if __name__ == "__main__":
    run_gui()
