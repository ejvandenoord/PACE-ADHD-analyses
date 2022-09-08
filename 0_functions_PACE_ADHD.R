#Function needed for collecting the files from RAMWAS
#Code written by Mandy Meijer, Karolina Aberg and Edwin van den Oord, 2020. Center for Biomarker Research and Precision Medicine (BPM), Virginia Commonwealth University
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




collect_results = function(param,analysis_label,collect_files) {

  # this RaMWAS comment adds the names of directories where RaMWAS will write results most relevant are param$dirpca and param$dirmwas
  param        = parameterPreprocess(param)

  results_dir = paste0(param$dirproject,"/",analysis_label,"/")
  dir.create(results_dir,showWarnings=F,recursive=T)
  
  collect_files = makefullpath(param,collect_files) # this function finds the exact location of the  files 

  file.copy( collect_files,results_dir) # this copies the files we need

  # this zips the files we need
  zip_file = paste0(results_dir,analysis_label,"_results.zip")
  zip::zip(zip_file,basename(collect_files), root=results_dir )
  message("Requested results archived in ",  zip_file )  
   
}


makefullpath = function(param,collect_files) {
  
  dirmwas_files = list.files(param$dirmwas,full.names=TRUE)
  dirpca_files  = list.files(param$dirpca,full.names=TRUE)
  dircov_files  = list.files(paste0(param$dirproject,"/covariates/"),full.names=TRUE)
  
  for (i in seq_along(collect_files)) { # i=1
    sel = grepl(collect_files[i],basename(dirpca_files) )
    if (any(sel)) collect_files[i] = dirpca_files[sel]
    
    sel = grepl(collect_files[i],basename(dirmwas_files) )
    if (any(sel)) collect_files[i] = dirmwas_files[sel]
    
    sel = grepl(collect_files[i],basename(dircov_files) )
    if (any(sel)) collect_files[i] = dircov_files[sel]
  }
  
  collect_files
}