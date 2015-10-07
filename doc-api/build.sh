#!/bin/bash

STATIC_NB_DIR=build/html/notebooks/
mkdir -p $STATIC_NB_DIR

for notebook in source/notebooks/*.ipynb
do
    echo $notebook $name
    name=${notebook%.ipynb}
    # Jupyter can convert directly to rst, yet there is an issue with the
    # titles <https://github.com/ipython/ipython/issues/8674>.
    # Therefore, we do the conversion in two steps.
    jupyter-nbconvert --to markdown --ExecutePreprocessor.enabled=True $notebook --output ${name}.md
    pandoc -i ${name}.md -o ${name}_.rst
    sed -Ee 's/:([^:]+):``([^`]+)``/:\1:`\2`/g' ${name}_.rst > ${name}.rst
    rm ${name}.md ${name}_.rst

    jupyter nbconvert --to html --ExecutePreprocessor.enabled=True $notebook --output $STATIC_NB_DIR/$(basename $name)_.html
    cp $notebook $STATIC_NB_DIR
done

make html
