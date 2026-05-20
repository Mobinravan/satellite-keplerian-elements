# Satellite Keplerian Elements Calculator

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Astrodynamics](https://img.shields.io/badge/Astrodynamics-Orbital_Mechanics-orange.svg)]()

Calculation of Keplerian and Hill's orbital elements for GEO and Molniya satellites from position and velocity vectors.

## Overview

This project calculates the complete set of orbital elements for two communication satellites:
- **INTELSAT-20** (Geostationary Earth Orbit - GEO)
- **MONILIYA3-50** (Molniya-type highly elliptical orbit)

Given position (r) and velocity (v) vectors in an inertial coordinate system, the program computes:
- 6 classical Keplerian elements
- True anomaly and Mean anomaly
- Hill's elements (r, ṙ, u, Ω, G, H) for the GEO satellite
- Orbital period verification using Kepler's 3rd law
- 2D orbit visualization in perifocal coordinates

## Problem Statement (Original)

> Given position and velocity vectors for INTELSAT-20 and MONILIYA3-50 satellites:
> 
> **INTELSAT-20:**
> - r = (16575.5879, -38775.5254, 5.3421) km
> - v = (2.8271, 1.2079, -0.0008) km/s
> 
> **MONILIYA3-50:**
> - r = (-495.4135, 8636.4268, -3136.0422) km
> - v = (-4.0917, 5.5089, 4.9534) km/s
> 
> Earth's gravitational constant (GM) = 398600.4418 km³/s² (WGS84)

## How to Run
```bash
python code/keplerian_elements.py
