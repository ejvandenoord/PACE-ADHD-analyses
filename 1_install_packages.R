#Installing required packages
#Code written by Mandy Meijer, Karolina Aberg and Edwin van den Oord, 2020. Center for Biomarker Research and Precision Medicine (BPM), Virginia Commonwealth University
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#This only needs to be done once

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "zip",
  "ramwas", 
  "minfi",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "wateRmelon",
  "readxl",
  "RPMM",
  "FlowSorted.Blood.450k",
  "FlowSorted.Blood.EPIC"),
  update = TRUE, ask = FALSE, quiet = TRUE)

