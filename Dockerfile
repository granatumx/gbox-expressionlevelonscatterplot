FROM granatumx/gbox-py-sdk:1.0.0

RUN pip install palettable

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
