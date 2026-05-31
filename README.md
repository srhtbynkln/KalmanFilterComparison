# Kalman Filter Variants for Railway Localization

This repository implements and compares multiple **Kalman Filter** architectures for the **Modena Train dataset**, comparing Kalman Filter algorithms on a railway localization problem by fusing GNSS and IMU.

## Implemented Algorithms
To evaluate different state estimation techniques, the following filters are implemented and compared in this repository:
- **KF** (Standard Kalman Filter)
- **Sage-Husa KF** (Sage Husa Adaptive Kalman Filter)
- **Strong Tracking KF** (Strong Tracking Adaptive Kalman Filter)
- **EKF** (Extended Kalman Filter)
- **Sage-Husa EKF** (Sage Husa Adaptive Extended Kalman Filter)


## Dataset

The **Modena Train dataset** used in this project is based on the following works:

**1. Amatetti et al., 2022**  
*Towards the Future Generation of Railway Localization and Signaling Exploiting sub-meter RTK GNSS*  
2022 IEEE Sensors Applications Symposium (SAS), pp. 1–6, IEEE.

**2. Mikhaylov et al., 2023**  
*Towards the Future Generation of Railway Localization Exploiting RTK and GNSS*  
IEEE Transactions on Instrumentation and Measurement, 2023, IEEE.

These publications provide the methodology and high-precision RTK GNSS data used to validate the EKF implementation.

## Features

- Fusion of GNSS and IMU
- Accurate estimation of train position and heading  

## Example Result
<img width="2089" height="1470" alt="graphics-1" src="https://github.com/user-attachments/assets/89e0af94-c24f-420b-af55-505c21bd35a9" />
<img width="1570" height="1106" alt="marmaray" src="https://github.com/user-attachments/assets/aa5ef389-237f-4df2-ba49-1e85eb5824c3" />
