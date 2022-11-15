> Trabajo práctico final | Introducción a la Bioinformática - UNQ
> Prof. Dra. Ana Julia Velez Rueda

# Anotación por homología de mirARNs

* [1 Objetivo](#1-objetivo)
* [2 Contexto](#2-contexto)
* [3 Requerimientos detallados](#3-requerimientos-detallados)
* [4 Tecnologías](#4-tecnologías)
* [5 Forma de entrega](#5-forma-de-entrega)

## [1 Objetivo](#1-objetivo)
El objetivo de este trabajo práctico es integrar los conceptos biológicos y bioinformáticos desarrollados durante la materia, junto con los conocimientos, prácticas y habilidades propias de la informática y la programación adquiridas en otras materias, a través de la construcción de un programa simple pero innovador que sirva de asistencia al proceso del análisis bioinformático y/o herramienta educativa para la enseñanza de la Biología. 

## [2 Contexto](#2-contexto)
Los micro-ARN son moléculas de ARN transcritas a partir de genes de ADN, pero no son traducidas a proteínas. Estos son en general ARN monocatenario, de una longitud de entre 21 y 25 nucleótidos, que tiene la capacidad de regular la expresión de otros genes mediante diversos procesos. Se expresan en una amplia variedad de organismos, desde plantas hasta mamíferos. Muchos micro-ARN están bien conservados entre especies,​ y muchos componentes de la maquinaria de los micro-ARN se han encontrado incluso en arqueas y bacterias, lo que revela que su origen es muy antiguo (Bartel 2004).
La secuencia de ADN que codifica para un gen de micro-ARN tiene una longitud que supera al tamaño final del propio micro-ARN e incluye la región micro-ARN y una región que es complementaria a la anterior, lo que permite su apareamiento. Esto conlleva que, durante la transcripción de esta secuencia de ADN, se forman regiones que tienen la capacidad de formar una horquilla y generar un ARN bicatenario primario largo conocido como pri-micro-ARN. Posteriormente, una enzima nuclear llamada DROSHA corta las bases de la horquilla, formando lo que se denomina pre-micro-ARN. Este pre-micro-ARN es transportado desde el núcleo al citoplasma por la exportina. Una vez que el pre-micro-ARN está en el citoplasma es fragmentado por la enzima DICER, que lo corta hasta su longitud final.
En la actualidad, existe una gran variedad de aplicaciones para el estudio y caracterización de miRNAs, que van desde el diseño de alimentos (Wani et al. 2020), hasta el tratamiento de enfermedades o regulación metabólica (Preethi and Sekar 2021; Chakraborty et al. 2017). 

## [3 Requerimientos detallados](#3-requerimientos-detallados)

El objetivo es desarrollar una biblioteca de Python que permita buscar de forma eficiente los miRNAs que regulan la expresión de un gen de referencia en un organismo dado:

1. El diseño del software deberá soportar tanto la entrada de una secuencia FASTA de un gen, como un Gene Id. Una corrida requerirá como input 1 secuencia FASTA (o Gene Id) de referencia y una especie target (nombre científico) sobre la cuál queremos anotar los miRNAs. 

2. Una vez iniciada la corrida, el software:
 Deberá disparar la búsqueda de genes homólogos utilizando [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) local,, con las bases de datos específicas para ácidos nucleicos. Sobre una lista ponderada de dichos hits, el software realizará la anotación de miRNAs para el gen homólogo del organismo target encontrado. La anotación de miRNA reguladores del gen se deberá realizar mediante el uso de la base de datos de alguna de las siguientes base de datos, configurable por la/el/le usuaria/o/e: [Rumir]( http://rumimir.sigenae.org/), [TarBase](https://www.hsls.pitt.edu/obrc/index.php?page=URL1237572545), [miRNEST](http://rhesus.amu.edu.pl/mirnest/copy/), [mirBase](https://www.mirbase.org/index.shtml)

> PLUS I: Posibilitar la anotación cruzando múltiples bases de datos, haciendo búsquedas en más de una o todas las antes mencionadas.

3. Deberá realizar la búsqueda de miRNA para el gen en cuestión, y sobre la lista de resultados obtenidos deberá disparar la búsqueda de miRNA homólogos a estos en en la especie target.

4. Sobre la base de los resultados obtenidos, el programa deberá permirmitir hacer un resumen de los resultados, integrando en una lista de resultados los miRNAs sin repeticiones que regulan el gen homólogo de la especie target y exportarlo en un formato conveniente.

> **REQUISITOS OBLIGATORIOS**: El programa deberá estar bien dodumentado para su uso con los comandos y parámetros que pueden ser utilizados. El proyecto en GitHub deberá contar con una wiki de documentación y archivos y tutoriales de  ejemplos de corrida. Construir una versión distribuída en paquete de PIP, que pueda ser instalada via Test Pypi con `pip install`.El software debe contar con tests unitarios.

## [4 Tecnologías](#4-tecnologías)
El programa podrá ser implementado utilizando Python. El programa deberá ser portable en los sistemas operativos Unix y Mac (Windows opcional). Algunas de las bibliotecas y herramientas recomendadas para la realización del trabajo son: 

* Python: argparse, Biopython
* PIP: https://pypi.org/ 

No son todas las herramientas deben ser usadas si o si y presentan características diferentes, además de que tampoco son las únicas opciones disponibles. Queda a criterio del equipo la elección de la biblioteca. 

## [5 Forma de entrega](#5-forma-de-entrega)
El trabajo práctico se realizará en equipos de hasta 4 integrantes, y deberá ser entregado de forma virtual el día 30  de Noviembre. El mismo deberá además estar en un repositorio público en Github, creado en la organización de la materia. Tod@s l@s integrantes del equipo deben tener commits con aportes significativos. Los trabajos serán presentados por los miembros del equipo en una exposición oral (con DEMO), la última clase de la cursada (2 de Diciembre). Los trabajos contarán con una calificación grupal por código y conceptual evaluando los contenidos de la materia. Esta nota contempla no sólo los aspectos de diseño del software y de apropiación de conceptos biológicos o de herramientas bioinformáticas, sino que también otros aspectos que se consideran importantes en la evaluación por proceso, como ser: entrega de versiones y consultas previas, indagación/investigación por fuera de lo propuesto por la docente, entregas en los tiempos pautados, trabajo colaborativo/en equipo. Cada estudiante tendrá además una calificación individual correspondiente a los aportes y su proceso personal.
