FROM conda/miniconda3
ADD . /code
RUN apt-get update -y
RUN apt-get install -y build-essential python-dev python-pmw libglew-dev freeglut3-dev libpng-dev libfreetype6-dev libxml2-dev libmsgpack-dev python-pyqt5.qtopengl libglm-dev libnetcdf-dev
RUN apt-get install -y git-all
RUN chmod 775 /code/docker_install.sh
RUN /code/docker_install.sh
RUN conda init bash
RUN cd /code; conda env create -f environment.yml 
