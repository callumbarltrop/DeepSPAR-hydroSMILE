# Clear the environment in R

rm(list=ls())

# Load in/install all of the required packages 

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

# Install python version 3.10.11. Make sure Git is installed prior to running the commands below. 

# For windows

py_version <- "3.10.11"
path_to_python <- reticulate::install_python(version=py_version)

# For linux (uncomment the lines below, change your username)

# module load release/24.04  GCC/12.3.0 R/4.3.2

# install_dir <- "/home/user/Python"
# path_to_python <- reticulate::install_python(
#   version = "3.10.11",
#   target = install_dir # <-- This is the key argument
# )

.rs.restartR() # Restart R 

# Re-load in packages 

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Install virtual environment 

reticulate::virtualenv_create(envname = 'DeepSPAR_env',
                              python=path_to_python,
                              version=py_version)

# Check we can use virtual environment 

reticulate::use_virtualenv("DeepSPAR_env", required = T)

.rs.restartR() # Restart R 

# Re-load in packages 

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Load virtual environment 

reticulate::use_virtualenv("DeepSPAR_env", required = T)

# Install tensorflow v2.11.0 into virtual environment 

tf_version="2.11.0" 
tensorflow::install_tensorflow(method="virtualenv", envname="DeepSPAR_env",
                               version=tf_version) 

.rs.restartR() # Restart R 

# Re-load in packages 

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Load virtual environment 

reticulate::use_virtualenv("DeepSPAR_env", required = T)

# Install keras into virtual environment 

keras::install_keras(method = c("virtualenv"), envname = "DeepSPAR_env",version=tf_version) #Install keras

# Install numpy into virtual environment 

reticulate::py_install("numpy==1.26.4", envname = "DeepSPAR_env", method = "virtualenv", pip = TRUE)

.rs.restartR() # Restart R 

# Re-load in packages 

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Load virtual environment 

reticulate::use_virtualenv("DeepSPAR_env", required = T)

# Check if Keras is available and working in the created virtual environment
# If you get an error here, or anywhere else in the script, us an LLM to debug with the errors you receive

keras::is_keras_available() # Should print 'TRUE'
