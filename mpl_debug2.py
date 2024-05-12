import streamlit as st
import pandas as pd
import plotly.graph_objects as go

def plot_graph(data):
    fig = go.Figure(data=[go.Scatter(x=data['x'], y=data['y'], z=data['z'], mode='markers')])
    st.plotly_chart(fig)

def main():
    points = [
        (30.00, 0.00, 0.00),
        (29.10, 6.76, 6.76),
        (23.11, 18.42, 18.42),
        (- 3.53, 37.73, 37.73),
        (- 106.35, 47.04, 47.04),
        (- 149.97, -30.32, -30.32),
        (- 142.39, -34.84, -34.84),
        (1.57, -35.27, -35.27),
        (9.43, -30.68, -30.68),
        (25.51, -14.96, -14.96),
        (28.33, -9.20, -9.20),
    ]
    st.title("Interactive Graph")

    if st.button("Add"):
        x = [int(points[point][0]) for point in points]
        y = [int(points[point][1]) for point in points]
        z = [int(points[point][2]) for point in points]
        points.extend(zip(x, y, z))

    if st.button("Plot"):
        data = pd.DataFrame(points, columns=["x", "y", "z"])
        plot_graph(data)


if __name__ == "__main__":
    main()
