ARG BASE_IMAGE=ubuntu:latest
FROM ${BASE_IMAGE}
RUN apt-get update
RUN apt-get install -y pip python3.12-venv
RUN apt-get install -y git 
RUN git clone --branch Removing-SetupTools https://github.com/ComboProblem/cutgeneratingfunctionology.git /user/cutgeneratingfunctionology
WORKDIR user/cutgeneratingfunctionology
ENV VIRTUAL_ENV=user/cutgeneratingfunctionology/venv-cgf
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install ".[passagemath]"
CMD ["bash"]
