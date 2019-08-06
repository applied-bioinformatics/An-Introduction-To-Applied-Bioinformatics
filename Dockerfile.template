# Distributed under the terms of the Modified BSD License.
#
# This image is expecting to be built alongside the *built* IAB docs,
# e.g. https://github.com/applied-bioinformatics/built-iab,
# and is designed for running on mybinder.org

FROM jupyter/minimal-notebook:58169ec3cfd3

LABEL maintainer="Greg Caporaso <gregcaporaso@gmail.com>"

USER root

ENV DISPLAY=:99

RUN apt-get update
RUN apt-get install -y xvfb x11-utils

COPY IAB-notebooks ${HOME}/IAB-notebooks/
COPY .jupyter/custom/custom.css ${HOME}/.jupyter/custom/custom.css
# `fix-permissions` ships with jupyter/minimal-notebook
RUN fix-permissions ${HOME}
RUN rm -rf work

USER ${NB_UID}

COPY environment.yml ${HOME}
RUN conda env update -n base -f environment.yml
RUN rm environment.yml

# This is almost identical to the `ENTRYPOINT` defined in jupyter/minimal-notebook,
# except we tack on a `xvfb-run` on the end, which ensures that ete3 has X.
ENTRYPOINT ["tini", "-g", "--", "xvfb-run"]
