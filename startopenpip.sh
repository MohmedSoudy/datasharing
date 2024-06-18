cd /home/aranjan/webapp/
docker-compose build && docker-compose up
touch docker.txt
docker ps >> docker.txt
input="docker.txt"
while IFS= read -r line
do
	b=${line:0:12}
	if [[ $line == *"mysql:8.0.0"* ]]; then
	echo "$b"; break
	fi
done < "$input"
rm docker.txt
c="/home/aranjan/webapp/db/dev10.0_huri.sql"
sudo docker cp  "$c" "$b":/
docker exec -it "$b" /bin/bash