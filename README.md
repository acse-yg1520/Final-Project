#  Bubble Analyser

Bubble Analyser is a Matlab Application presented by Mesa et al. (2022) for processing and quantifying bubble images, normally taken for assessing froth flotation equipment. The one dimensional Discrete Fourier Transform (Vinnett et.al, 2018), circles detection method (Wang et al., 2022) combined with Watershed segmentation and partial Watershed are implemented to estimate the Sauter mean diameter of bubbles are implemented in this software. For more details about installation and usage, you can refer to `Bubble Analyser Manual.pdf`.

## Installation

You can simply download a copy from this repository as a .zip (or equivalent) file. Unzip the file and keep the unpacked folder wherever you prefer—We recommend to put it in the default Matlab folder. Once you have decided where to keep it, make sure you add this location to your Matlab path. For help on doing this, please refer to this [link](https://www.mathworks.com/help/matlab/ref/addpath.html).

## Usage

If you have added the ‘Bubble Analyser’ folder to your Matlab path, you can run the APP by:
```>> BubbleAnalyserApp()```

Bubble Analyser was developed for Matlab 2020, or newer. However, it might run with older versions of Matlab. Bubble Analyser uses the following Toolboxes: Image Processing and Signal Processing.

You can find the sample bubble images inside the folder `Sample_photos`.

## Reference
Mesa, Diego, Quintanilla, Paulina, and Reyes, Francisco (2022). Bubble Analyser — An open-source software for bubble size measurement using image analysis. Minerals Engineering, Volume 180, 2022, 107497, ISSN 0892-6875 https://doi.org/10.1016/j.mineng.2022.107497.

L.Vinnett, J. Sovechles, C.O. Gomez, K.E. Waters, An image analysis approach to determine average bubble sizes using one-dimensional Fourier analysis, Minerals Engineering, Volume 126, 2018, Pages 160-166, ISSN 0892-6875, https://doi.org/10.1016/j.mineng.2018.06.030.

Wang, J.; Forbes, G.; Forbes, E. Frother Characterization Using a Novel Bubble Size Measurement Technique. Applied. Sci. 2022, 12, 750. https://doi.org/10.3390/app12020750

## License
Bubble Analyser is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 3 only of the License.

Bubble Analyser is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Bubble Analyser.  If not, see <https://www.gnu.org/licenses/>.
