<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All experiments for uploaded samples
        </h3>
        <table align="center" class="colored">
            %if len(experiments) == 0:
                <tr>
                    <td colspan="2">
                        There are no experiments for uploaded samples
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Seq Facility</td>
                    <td>Array Version</td>
                    <td>Result Folder Name</td>
                    <td>Plate Barcode</td>
                    <td>Samples of this Experiment</td>
                </tr>
                <% ctr = 0 %>
                %for experiment in experiments:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${experiment[1]}</td>
                        <td>${experiment[2]}</td>
                        <td>${experiment[3]}</td>
                        <td>${experiment[4]}</td>
                        <td><a href="${h.url_for( controller='samples', action='of_experiment', sort_id='default', order='default', experiment_id=experiment[0], seq_facility=experiment[1], array_version=experiment[2], result_folder_name=experiment[3], plate_barcode=experiment[4] )}">Samples</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

