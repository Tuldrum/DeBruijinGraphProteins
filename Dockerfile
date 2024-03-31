    # Usar una imagen base de Python
    FROM python

    WORKDIR /app

    # Copiar el archivo de requerimientos e instalar las dependencias
    COPY requirements.txt requirements.txt
    RUN pip install -r requirements.txt

    # Copiar los archivos de la aplicación a la imagen
    COPY . .

    # Especificar el comando a ejecutar cuando se inicie el contenedor
    # CMD ["python", "/proccesor_app/Proccesor.py"]