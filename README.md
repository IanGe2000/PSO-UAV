# PSO-UAV

Code of my graduate thesis for bachelor degree at BIT.  
UAV swarm trajectory planner using modified particle swarm optimization (PSO).

## In this repo

|directory|description|
|-|-|
|```/Octave_simulation```|A simulation of the algorithm on GNU Octave.|
|```/CXX_simulation```|A C++ implementation with Eigen-based data structure of the algorithm.|
|```/PSO-UAV_Input```|A Qt-based GUI for parameter input and algorithm presentation.|

## Dependencies

### Octave

Install Octave for GNU/Linux,

```bash
sudo apt install octave
```

and start octave by

```bash
octave
```

### CMake

```bash
sudo apt install cmake
sudo apt install build-essential
```

### Eigen  

Get Eigen's source code, and I moved what I used of Eigen to ```/usr/local/include/```.  
If you have installed Eigen in other fashion before, or want to put the source code elsewhere, you will need to tell the compiler how to find them when compiling.

```bash
git clone https://gitlab.com/libeigen/eigen.git
sudo cp eigen/Eigen /usr/local/include/ -r
sudo cp eigen/unsupported/Eigen /usr/local/include/Eigen/unsupported -r
```

### Qt

First, download Qt 5.9.9 installer and run

```bash
wget https://download.qt.io/archive/qt/5.9/5.9.9/qt-opensource-linux-x64-5.9.9.run
chmod +x qt-opensource-linux-x64-5.9.9.run
./qt-opensource-linux-x64-5.9.9.run 
```

In the installation process, if you encountered

```text
Warning: Network error: [ QNetworkReply::NetworkError(AuthenticationRequiredError) ] "Authentication failed."
```

try Settings(bottom-left of the installation window)->no proxy.  

Next, at "Select Components" page, this repo only need ```Qt 5.9.9/Sources``` and ```Qt 5.9.9/Qt Charts```.  
Finally, for Qt 5.9.9 an additional package is needed for compiling.

```bash
sudo apt install libqt5charts5-dev
```

## Compilation guide

```bash
cd PSO-UAV
mkdir build
cd build
cmake ..
cmake --build .
```

Now, you will have the C++ implementation of the algorithm built under ```build/CXX_simulation``` as ```CXX_simulation```, and the Qt-based GUI built under ```build/PSO-UAV_Input``` as ```PSO-UAV_Input```.  

For GUI presentation please follow

```bash
cp CXX_simulation/CXX_simulation PSO-UAV_Input/
./PSO-UAV_Input/PSO-UAV_Input
```
