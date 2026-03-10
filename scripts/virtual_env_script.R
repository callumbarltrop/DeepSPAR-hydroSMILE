packages = c("keras","tensorflow","reticulate")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

py_version <- "3.10.11"
path_to_python <- reticulate::install_python(version=py_version)

#For linux 

# module load release/24.04  GCC/12.3.0 R/4.3.2

# install_dir <- "/home/camu798f/Python"
# path_to_python <- reticulate::install_python(
#   version = "3.10.11",
#   target = install_dir # <-- This is the key argument
# )



reticulate::virtualenv_create(envname = 'DeepSPAR_env',
                              python=path_to_python,
                              version=py_version)

reticulate::use_virtualenv("DeepSPAR_env", required = T)

tf_version="2.11.0" 
reticulate::use_virtualenv("DeepSPAR_env", required = T)
tensorflow::install_tensorflow(method="virtualenv", envname="DeepSPAR_env",
                               version=tf_version) #Install version of tensorflow in virtual environment

keras::install_keras(method = c("virtualenv"), envname = "DeepSPAR_env",version=tf_version) #Install keras

reticulate::py_install("numpy==1.26.4", envname = "DeepSPAR_env", method = "virtualenv", pip = TRUE)

reticulate::use_virtualenv("DeepSPAR_env", required = T)
keras::is_keras_available() #Check if keras is available
