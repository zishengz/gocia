

![Picture1](./gocia_logo.png)

# G O C I A

**G**lobal **O**ptimizer for **C**lusters, **I**nterfaces, and **A**dsorbates

**GOCIA** es un kit de herramientas de optimización global y módulos de Python especializados en el muestreo de clusters soportados, interfaces reestructuradas y configuraciones de adsorbatos.

Copyright © 2020 Zisheng Zhang

Por favor, CITAR el siguiente artículo si utiliza cualquier parte de este repositorio:

> Zhang, Z.; Alexandrova, A. N., GOCIA: a grand canonical global optimizer for clusters, interfaces, and adsorbates. Phys. Chem. Chem. Phys., 2025,27, 696-706. doi:[10.1039/D4CP03801K](https://doi.org/10.1039/D4CP03801K). -> [PDF DOWNLOAD](https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/66fe375651558a15efeab3a6/original/gocia-grand-canonical-global-optimizer-for-clusters-interfaces-and-adsorbates.pdf)

[TOC]

## Requisitos

- Python 3.6 o posterior
- ASE y sus dependencias
- Natsort y LATEX (generación de informes en PDF)

## Instalación

### Entorno de Python 
Primero, instale su propio entorno de Python, ya que las HPC (supercomputadoras) generalmente no otorgan a los usuarios regulares permisos de escritura en la ruta de Python. Para ahorrar espacio en disco, se recomienda instalar ```Miniconda```:

```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

Acepte todas las opciones por defecto y ejecute ```source ~/.bashrc``` para activar el entorno de conda. Puede ejecutar ```python``` en la terminal para verificar la versión de Python que está utilizando.


### Instalar `GOCIA`

Si su máquina tiene `Git` instalado, simplemente clone el repositorio en su directorio local mediante:

```bash
git clone https://github.com/zishengz/gocia.git
```

O, también puede descargar y descomprimir el código fuente:

```bash
wget https://github.com/zishengz/echo/archive/refs/heads/main.zip
unzip main.zip
rm main.zip
mv main gocia
```

Una vez obtenido el repositorio `gocia`, agréguelo a su `PYTHONPATH` mediante:

```bash
export PYTHONPATH=$PYTHONPATH:`pwd`/gocia/
```

Recuerde agregar esta línea de exportación a su `~/.bashrc` o al script de envío, de modo que el paquete `GOCIA` sea accesible por Python al ejecutar el trabajo.

Debe usar la ruta absoluta (puede verificarla ejecutando `pwd` en la terminal Bash) para este propósito.

Después de esto, ejecute la siguiente línea para probar:

```bash
python -c 'import gocia'
```
Si no aparece ningún error, ¡GOCIA debería haberse importado a su ruta!

### Actualizar `GOCIA`

Si instaló mediante Git, actualice extrayendo los cambios desde la rama main en el directorio `gocia`:

```bash
cd xxx/gocia
git pull
```

De lo contrario, debe eliminar manualmente el directorio antiguo `gocia` y luego volver a descargar y descomprimir.


## Tutorial

Asumimos el uso de VASP para la optimización local, a menos que se especifique lo contrario.

```HPC``` representa el planificador de trabajos en el clúster que utiliza:

- slurm: CORI
- sge: Hoffman2
- pbs: DoD machines

### Optimización local en 3 pasos

Los archivos necesarios incluyen:

- INCAR-1, INCAR-2, INCAR-3
  para cálculos DFT de baja, media y alta precisión.
- KPOINTS
- init-worker.py
  El trabajo del trabajador que ejecuta optimizaciones locales en 3 pasos y verifica la conectividad poco razonable.
- HPC-vasp-init.sh
  El script de shell para enviar trabajos del trabajador.
  RECUERDE reemplazar la ruta de ```.bashrc``` por la suya.
- input.py
  Un archivo de datos que contiene la ruta del pseudopotencial y el comando de VASP.
- substrate.vasp (opcional)
  Para referencia durante la verificación de la geometría, si se proporciona zLim.
  

Procedimiento:

1. Reemplace la ruta de ```.bashrc``` en ```HPC-vasp-init.sh``` por la suya.
2. Coloque la ruta a sus pseudopotenciales y el comando de VASP en ```input.py```
3. Conceda permisos de ejecución al script de envío mediante ```chmod +x HPC-vasp-init.sh```
4. Envíe el trabajo mediante ```./HPC-vasp-init.sh xxx.vasp```

### Población inicial: Generación estructural

Los archivos necesarios incluyen:

- substrate.vasp
  El archivo de estructura en formato VASP que contiene la losa del sustrato, con restricciones.
- xxxSample.py
  Elija el método de muestreo estructural que mejor se adapte a su sistema.

```bash
python xxxSample.py substrate.vasp
```

### Población inicial: Optimización local

Los archivos necesarios incluyen:

- INCAR-1, INCAR-2, INCAR-3
- KPOINTS
- init-worker.py
- HPC-vasp-init.sh
- input.py
- db2vasp.py
  Script para convertir archivos de base de datos ase a archivos en formato VASP con nombres sistemáticos.
- collectVASP.py
  Escribe los resultados de VASP en un archivo de base de datos ase y filtra los duplicados.

Procedimientos (1-3 son los mismos que en la sección de optimización en 3 pasos):

1. Reemplace la ruta de ```.bashrc``` en ```HPC-vasp-init.sh``` por la suya.
2. Coloque la ruta a sus pseudopotenciales y el comando de VASP en ```input.py```
3. Conceda permisos de ejecución al script de envío mediante ```chmod +x HPC-vasp-init.sh```
4. Convierta el archivo .db a formato VASP mediante ```python db2vasp.py xxx.db```
5. Envíe en lote por: ```for i in s0*vasp; do ./HPC-vasp-init.sh $i; done```
6. Después de que finalicen todos los trabajos, recoja los resultados mediante ```python collectVASP.py```

### Muestreo GCGA

Los archivos necesarios incluyen:

- substrate.vasp
- INCAR-1, INCAR-2, INCAR-3
- KPOINTS
- ga-HPC.py
  El trabajo principal que se ejecuta localmente (en el nodo de inicio de sesión si está permitido) y controla el envío de trabajos.
- ga-worker.py
  El trabajo del trabajador que ejecuta optimizaciones locales y actualiza la población en los nodos de cálculo.
- HPC-vasp.sh
  El script de shell para enviar trabajos del trabajador GCGA.
- input.py
  Un archivo de datos que contiene la información necesaria para el muestreo GCGA.
- gcga.db
  El archivo de base de datos ase que contiene la población inicial, obtenido del paso anterior.

Procedimientos:

1. Reemplace la ruta de ```.bashrc``` en ```HPC-vasp.sh``` por la suya.
2. Coloque la ruta a sus pseudopotenciales, el comando de VASP, los potenciales químicos y otros parámetros GCGA en ```input.py```
3. Conceda permisos de ejecución al script de envío mediante ```chmod +x HPC-vasp.sh```
4. Copie el archivo de base de datos de [población inicial: optimización local] a ```gcga.db``` 
5. Ejecute el maestro GCGA en el nodo de inicio de sesión mediante ```nohup python -u ga-HPC.py &```
6. Si desea detener el GCGA, ejecute ```touch STOP```.



Otras partes están en construcción...
