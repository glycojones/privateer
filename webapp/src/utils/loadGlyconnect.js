async function get_glyconnect_id(glytoucan_id) {
    let url = "https://glyconnect.expasy.org/api/structures/search/glytoucan"

    console.log(JSON.stringify({"glytoucan_id": glytoucan_id}))

    fetch(url, {
        method: "POST",
        body: JSON.stringify({"glytoucan_id": glytoucan_id}),
        headers: {
            "Content-type": "application/json; charset=UTF-8"
        }
    }).then((response) => {
        const contentType = response.headers.get("content-type")
        if (contentType && contentType.indexOf("application/json") !== -1) {
            return response.json()
        } else if (response.status == 404) {
            return Promise.reject("404 Error")
        }
    }).catch((error) => {
        return Promise.reject(error)
        // console.log(error)
    })
}

export default async function load_glyconnect(table_data) {

    let promises = [];

    for (var i = 0; i < table_data.length; i++) {
        let item = table_data[i]
        promises.push(get_glyconnect_id(item.glytoucan_id))
    }

    try {
        const data = await Promise.all(promises)
    } catch {

    }
    data.forEach((data, index) => {
        // console.log(data)
        // table_data[index].glytoucan_id = data[0].id;
    })
    return table_data;
}
