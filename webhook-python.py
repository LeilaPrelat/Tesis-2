#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

mandar msj a discord cuando un script de python termina

"""

from discord import Webhook, RequestsWebhookAdapter

url = 'https://discord.com/api/webhooks/798651656202878986/ViaIJQ9vQOCZa2U2NbIFYLvOTrNF4vID5GCW5_iB9Ozu1VA0edqH4B7a0Uot2v_syYn-'

webhook = Webhook.from_url(url, adapter=RequestsWebhookAdapter())
webhook.send('python termino')
