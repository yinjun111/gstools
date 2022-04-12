#synchronize scripts in the dev codes folder with the running folder
rsync -a -d --delete /home/jyin/Pipeline/gstools/ /apps/gstools/

/usr/local/bin/aws s3 sync /home/jyin/Pipeline/gstools/ s3://ferring-omictools/gstools/ --delete

#make everything execuable +x
chmod 755 -R /apps/gstools/
