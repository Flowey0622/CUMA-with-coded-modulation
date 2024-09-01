# CUMA-with-coded-modulation

## Description
Due to the recent growth in mobile data demand, the sixth generation (6G) requires higher throughput and coverage than the fifth generation (5G). To achieve this goal, we propose a new multiple access technique named Compact Ultra Massive Antenna Array (CUMA) to realise massive connectivity. In this simulation, the antennas are spatially tightly spaced at the transmitter side and a two-dimensional Fluid Antenna System (2D-FAS) is used at each user equipment (UE), which requires a large number of ports to be selected for reception during the channel transformation process, which are aggregated to be the signals that need to be decoded. Specifically, the UE only needs two radio-frequency (RF) chains to constructively combine the in-phase and quadrature components of the desired signal, while the interference signals are randomly accumulated. To enhance the reliability and efficiency of original CUMA systems, Turbo codes and LDPC codes are added, which improve error correction by introducing redundancy, allowing the system to recover error data during transmission. Our results will demonstrated based on the BER which is the main metric to evaluate the CUMA system's performance. 

### Comparisons  
1. modulation schemes: BPSK, QPSK, 16-QAM, 64-QAM and 256-QAM
2. different parameters: the number of multipath L,  the rice factor K in the Rician channel, the number of selected ports, the total number of ports N, the SNR, and the number of UEs

## Project structure
- Original CUMA System
- CUMA with Turbo codes
- CUMA with LDPC codes
- README.md
```
- `Original CUMA System/`: The base version of CUMA system under different modulations and compare their BER with different parameters.
- `CUMA with Turbo codes/`: The CUMA system with Turbo codes contains a performance comparison under different modulation schemes and parameters. The post-run data has been saved to a /data for further analysis.
- `CUMA with LDPC codes/`: The CUMA system with LDPC codes contains a performance comparison under different modulation schemes and parameters. The post-run data has been saved to a /data for further analysis.
- `README.md`:The README file for the project.

## Program run instruction
1. Clone the repository
2. Open the matlab and navigate to the code directory
3. Open each folder to run the corresponding part
4. The instructions for each part are detailed in the readme file within each folder
5. Execute the provided scripts to generate data and compare BER performance
6. Analyze the data: Review the generated data in the /data folder to compare the performance under different conditions
