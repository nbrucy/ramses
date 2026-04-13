if [ ! -e "ramses.data" ]; then
    echo "ramses.data does not exist, dowloading it..."
    wget --timeout=10 --tries=3 --no-check-certificate https://farou.ynh.fr/nextcloud/s/b5j6KMdcno6go2Z/download
    mv download ramses.data


    echo "ramses.data downloaded."
else
    echo "ramses.data already exists."
fi
