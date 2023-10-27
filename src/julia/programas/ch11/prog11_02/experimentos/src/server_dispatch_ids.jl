using HTTP.WebSockets

global id_valor
id_valor = 1

server = WebSockets.listen!("127.0.0.1", 8081) do ws

        for msg in ws
            send(ws, id_valor)
            global id_valor += 1
        end
    end

close(server)