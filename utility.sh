
check_mkdir()
{
if [ -d $1 ]
then
    echo "$1 exists."
else
    echo "$1 does not exits."
    mkdir -p $1
fi
}

check_error()
{
if [ $1 -ne 0 ]; then
    echo "FATAL ERROR: pipeline script"
    echo "ERROR CODE: $1"
    exit $1
fi
}

check_file_exists()
{
if [ -f $1 ]; then
  echo "$1 exists."
  return
fi
echo "$1 does not exists."
exit 1
}
