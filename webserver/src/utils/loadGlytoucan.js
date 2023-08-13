async function get_glytoucan_id(wurcs) {
    let url = "https://api.glycosmos.org/sparqlist/wurcs2gtcids?wurcs=" + encodeURIComponent(wurcs)

    return fetch(url, {
        method: "GET" // default, so we can ignore
    })
        .then((response) => response.json())
        .catch((error) => console.log(error))
}

export default async function load_glytoucan(table_data) {

    let promises = [];

    for (var i = 0; i < table_data.length; i++) {
        let item = table_data[i]
        // const id = get_glyconnect_id(item.wurcs)
        promises.push(get_glytoucan_id(item.wurcs))
    }

    const data = await Promise.all(promises)

    data.forEach((data, index) => {
        table_data[index].glytoucan_id = data[0].id;
    })
    return table_data;
}